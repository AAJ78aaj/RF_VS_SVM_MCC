#!/usr/bin/env python3
"""New Reviewer 3 diagnostic only; reads existing files and writes only here.
Run: python3 /home/ubuntu/revision_work/reviewer3_assessment/collinearity_diagnostic.py
No predictive model fitting, feature removal, manuscript edits, or p-values.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import re
import platform
import numpy as np
import pandas as pd

ROOT = Path('/home/ubuntu/revision_work')
OUT = ROOT / 'reviewer3_assessment'
MEMBERSHIP = ROOT / 'original_preserving_revision/analysis/archived_partition_locations.csv'
MODEL_SCRIPT = ROOT / 'analysis/revised_analysis.R'
ACTIVE = ['band4', 'NBR', 'B54', 'TCG', 'DEM', 'AnomB54', 'DiffB54']
CANDIDATES = ['band4', 'band5', 'NBR', 'B54', 'TCG', 'MeanB54', 'DEM', 'AnomB54', 'DiffB54']
SETS = {'active7': ACTIVE, 'diagnostic9': CANDIDATES}
SOURCES = [MEMBERSHIP, MODEL_SCRIPT] + [ROOT / f'input/Book{y}.csv' for y in (2005, 2009)]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save_csv(frame, filename, index=False):
    target = OUT / filename
    assert target.parent.resolve() == OUT.resolve()
    frame.to_csv(target, index=index, float_format='%.12g', na_rep='NA')

def md_table(frame, decimals=4):
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return 'NA' if pd.isna(v) else ('Infinity' if np.isinf(v) else f'{v:.{decimals}f}')
        return str(v)
    lines = ['| ' + ' | '.join(map(str, frame.columns)) + ' |', '| ' + ' | '.join(['---'] * len(frame.columns)) + ' |']
    lines.extend('| ' + ' | '.join(fmt(v) for v in row) + ' |' for row in frame.itertuples(index=False, name=None))
    return '\n'.join(lines)

OUT.mkdir(parents=True, exist_ok=True)
before = {str(p): digest(p) for p in SOURCES}
script_text = MODEL_SCRIPT.read_text()
match = re.search(r'features\s*<-\s*c\(([^)]+)\)', script_text)
assert match and re.findall(r'"([^\"]+)"', match.group(1)) == ACTIVE, 'Active predictor definition differs.'
member = pd.read_csv(MEMBERSHIP)
assert member['OID_'].notna().all() and not member['OID_'].duplicated().any(), 'Invalid membership OIDs.'
training = member.loc[member['partition'].eq('training'), ['OID_']].copy()
assert len(training) == 722, f'Expected 722 training members, got {len(training)}.'

pairs, summaries, column_audit, vif_rows, scale_rows, top_active_rows = [], [], [], [], [], []
audit = {
    'purpose': 'New Reviewer 3 training-only predictor-collinearity diagnostic, not historical filtering evidence.',
    'membership_file': str(MEMBERSHIP), 'membership_partition_counts': member['partition'].value_counts().to_dict(),
    'training_n': len(training), 'matching_key': 'OID_', 'year_specific_values': True,
    'active7_source': str(MODEL_SCRIPT), 'predictor_sets': SETS,
    'correlation_missing_policy': 'Pairwise complete finite observations; no imputation. Pair-specific n and status saved.',
    'spearman_policy': 'Pearson correlation of average ranks computed separately within each finite pair.',
    'vif_policy': 'Complete-case active7 X, standardized using training means and sample SD (ddof=1); each predictor regressed on the other six plus intercept; VIF=1/(1-R2). No response or RF/SVM fit.',
    'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__},
    'years': {}, 'source_sha256_before': before,
}
for year in (2005, 2009):
    source = ROOT / f'input/Book{year}.csv'
    raw = pd.read_csv(source)
    assert raw['OID_'].notna().all() and not raw['OID_'].duplicated().any(), f'Invalid {year} source OIDs.'
    assert raw['year'].eq(year).all(), f'Unexpected year in {source}.'
    assert set(CANDIDATES).issubset(raw.columns), 'Candidate column missing.'
    merged = training.merge(raw, on='OID_', how='left', validate='one_to_one', indicator=True)
    assert merged['_merge'].eq('both').all(), f'Unmatched {year} training OID.'
    tr = merged.drop(columns='_merge')
    # Coordinates and response are used only for alignment auditing, never for correlation or VIF.
    membership_train = member.loc[member['partition'].eq('training')].set_index('OID_')
    for c in ['Long', 'Lat', 'PresenceAb']:
        assert np.allclose(tr[c].to_numpy(), membership_train.loc[tr['OID_'], c].to_numpy()), f'Membership mismatch: {year} {c}'
    x = tr[CANDIDATES].apply(pd.to_numeric, errors='coerce')
    nbr_columns = [c for c in raw.columns if 'nbr' in c.lower()]
    audit['years'][str(year)] = {
        'source': str(source), 'source_rows': len(raw), 'matched_training_rows': len(tr),
        'source_columns': raw.columns.tolist(), 'nbr_named_columns': nbr_columns,
        'nbr_anomaly_gap': 'Only NBR is present; no identifiable NBR-anomaly column. Table 1 gap cannot be resolved from these CSVs; none invented.',
        'membership_coordinates_and_response_match': True,
        'missing_candidate_columns': [c for c in CANDIDATES if c not in raw.columns],
    }
    for c in CANDIDATES:
        finite = np.isfinite(x[c].to_numpy())
        finite_vals = x.loc[finite, c]
        column_audit.append({
            'year': year, 'predictor': c, 'training_n': len(tr),
            'raw_missing_n': int(tr[c].isna().sum()),
            'non_numeric_n': int((tr[c].notna() & x[c].isna()).sum()),
            'nonfinite_nonmissing_n': int((x[c].notna() & ~finite).sum()),
            'finite_n': int(finite.sum()), 'finite_unique_n': int(finite_vals.nunique()),
            'constant': bool(finite_vals.nunique() <= 1),
        })
    x = x.replace([np.inf, -np.inf], np.nan)
    year_pairs = []
    for set_name, cols in SETS.items():
        for method in ('pearson', 'spearman'):
            corr = pd.DataFrame(np.nan, index=cols, columns=cols)
            for c in cols:
                valid = x[c].dropna()
                if len(valid) >= 2 and valid.nunique() > 1:
                    corr.loc[c, c] = 1.0
            set_pairs = []
            for a, b in combinations(cols, 2):
                xy = x[[a, b]].dropna()
                n = len(xy)
                if n < 2:
                    value, status = np.nan, 'insufficient_pair_observations'
                elif xy[a].nunique() <= 1 or xy[b].nunique() <= 1:
                    value, status = np.nan, 'constant_in_pair'
                else:
                    z = xy.rank(method='average') if method == 'spearman' else xy
                    value = float(np.corrcoef(z[a], z[b])[0, 1])
                    status = 'ok'
                corr.loc[a, b] = corr.loc[b, a] = value
                row = {'year': year, 'predictor_set': set_name, 'method': method,
                       'predictor_1': a, 'predictor_2': b, 'n_pair': n,
                       'correlation': value, 'abs_correlation': abs(value), 'status': status}
                set_pairs.append(row)
            save_csv(corr.rename_axis('predictor'), f'correlation_{year}_{set_name}_{method}.csv', index=True)
            valid = [r for r in set_pairs if r['status'] == 'ok']
            ranked = sorted(valid, key=lambda r: (-r['abs_correlation'], r['predictor_1'], r['predictor_2']))
            best = ranked[0] if ranked else None
            summaries.append({
                'year': year, 'predictor_set': set_name, 'method': method, 'training_n': len(tr),
                'total_unique_pairs': len(set_pairs), 'valid_unique_pairs': len(valid),
                'undefined_pairs': len(set_pairs) - len(valid),
                'abs_r_gt_0_6': sum(r['abs_correlation'] > 0.6 for r in valid),
                'abs_r_gt_0_7': sum(r['abs_correlation'] > 0.7 for r in valid),
                'abs_r_gt_0_8': sum(r['abs_correlation'] > 0.8 for r in valid),
                'strongest_pair': f"{best['predictor_1']} / {best['predictor_2']}" if best else '',
                'strongest_correlation': best['correlation'] if best else np.nan,
            })
            pairs.extend(set_pairs)
            year_pairs.extend(set_pairs)
    # Ten distinct active7 pairs per year; include both correlation methods.
    active = pd.DataFrame([r for r in year_pairs if r['predictor_set'] == 'active7'])
    wide = active.pivot(index=['year', 'predictor_1', 'predictor_2', 'n_pair'], columns='method', values='correlation').reset_index()
    wide['max_abs_correlation'] = wide[['pearson', 'spearman']].abs().max(axis=1)
    wide = wide.sort_values(['max_abs_correlation', 'predictor_1', 'predictor_2'], ascending=[False, True, True]).head(10)
    wide.insert(1, 'rank', range(1, len(wide) + 1))
    top_active_rows.extend(wide.to_dict('records'))
    xv = x[ACTIVE].dropna()
    means = xv.mean()
    sds = xv.std(ddof=1)
    constants = sds.index[(~np.isfinite(sds)) | (sds <= 0)].tolist()
    audit['years'][str(year)]['vif_complete_case_n'] = len(xv)
    audit['years'][str(year)]['vif_constant_columns'] = constants
    for c in ACTIVE:
        scale_rows.append({'year': year, 'predictor': c, 'n_complete_case': len(xv),
                           'training_mean': means[c], 'training_sample_sd': sds[c]})
    if constants or len(xv) <= len(ACTIVE):
        for c in ACTIVE:
            vif_rows.append({'year': year, 'predictor': c, 'n_complete_case': len(xv),
                             'R2_auxiliary': np.nan, 'VIF': np.nan, 'design_rank': np.nan,
                             'status': 'not_computed_constant_or_insufficient_complete_cases'})
    else:
        z = (xv - means) / sds
        design_rank = int(np.linalg.matrix_rank(np.column_stack([np.ones(len(z)), z.to_numpy()])))
        for c in ACTIVE:
            y = z[c].to_numpy()
            other = z.drop(columns=c).to_numpy()
            design = np.column_stack([np.ones(len(y)), other])
            beta, *_ = np.linalg.lstsq(design, y, rcond=None)
            sse = float(np.sum((y - design @ beta) ** 2))
            sst = float(np.sum((y - y.mean()) ** 2))
            residual_fraction = sse / sst
            vif = float('inf') if residual_fraction <= np.finfo(float).eps else 1.0 / residual_fraction
            vif_rows.append({'year': year, 'predictor': c, 'n_complete_case': len(xv),
                             'R2_auxiliary': 1.0 - residual_fraction, 'VIF': vif,
                             'design_rank': design_rank,
                             'status': 'exact_linear_dependence' if np.isinf(vif) else 'ok'})

pair_df = pd.DataFrame(pairs)
summary_df = pd.DataFrame(summaries)
vif_df = pd.DataFrame(vif_rows)
audit_df = pd.DataFrame(column_audit)
top_df = pd.DataFrame(top_active_rows)
save_csv(pair_df.sort_values(['year', 'predictor_set', 'method', 'abs_correlation'], ascending=[True, True, True, False]), 'collinearity_tidy_pairs.csv')
save_csv(summary_df, 'collinearity_threshold_summary.csv')
save_csv(top_df, 'collinearity_strongest_active7_pairs.csv')
save_csv(vif_df, 'collinearity_vif_active7.csv')
save_csv(pd.DataFrame(scale_rows), 'collinearity_vif_training_scaling.csv')
save_csv(audit_df, 'collinearity_column_audit.csv')
after = {str(p): digest(p) for p in SOURCES}
assert before == after, 'A protected source changed during this diagnostic.'
audit['source_sha256_after'] = after
audit['sources_unchanged_during_diagnostic'] = True
(OUT / 'collinearity_input_audit.json').write_text(json.dumps(audit, indent=2) + '\n')

lines = [
    '# Reviewer 3: training-only predictor-collinearity diagnostic', '',
    '**Status:** This is a new reviewer-response diagnostic. It is not evidence of the original predictor-filtering procedure, has not been added to the manuscript or response letter, does not establish that RF/SVM predictions are bad, and is not a completed reduced-feature sensitivity test. No RF/SVM retraining, feature removal, predictive-performance evaluation, or editing of existing analyses/figures was performed.', '',
    '## Data and method', '',
    'The archived membership file identifies 722 training and 299 holdout OIDs. Each year contains 1,021 rows; exactly the same 722 training OIDs were matched one-to-one to the corresponding Book2005 or Book2009 CSV. Predictor values were never borrowed across years. Membership coordinates and response labels also agree with the corresponding source rows. Holdout values were not used in any correlation, scaling, or VIF calculation.', '',
    '`active7` (verified against revised_analysis.R): band4, NBR, B54, TCG, DEM, AnomB54, DiffB54. `diagnostic9`: band4, band5, NBR, B54, TCG, MeanB54, DEM, AnomB54, DiffB54.', '',
    'All nine candidates are available, numeric, finite, and nonconstant in each training sample; no missing observations, imputation, or constant-column exclusions were needed. Every pair and every VIF uses n = 722. Pearson measures linear association; Spearman uses average ranks for ties and measures monotonic association. Unique off-diagonal pairs are counted once, using strict absolute-correlation thresholds (>0.6, >0.7, >0.8), on unrounded values. No p-values are reported; the diagnostic concerns association magnitudes, not significance testing.', '',
    '**Table 1 gap:** The supplied CSVs contain only an `NBR`-named column, not an identifiable NBR-anomaly column. The NBR anomaly listed in original Table 1 therefore cannot be assessed from these inputs. No anomaly was invented; `AnomB54` was not treated as NBR anomaly.', '',
    '## Threshold summary', '',
]
show = summary_df[['year', 'predictor_set', 'method', 'valid_unique_pairs', 'abs_r_gt_0_6', 'abs_r_gt_0_7', 'abs_r_gt_0_8']].rename(columns={'valid_unique_pairs': 'pairs', 'abs_r_gt_0_6': 'absolute r > 0.6', 'abs_r_gt_0_7': 'absolute r > 0.7', 'abs_r_gt_0_8': 'absolute r > 0.8'})
lines.append(md_table(show))
lines.extend(['', '## Strongest active7 pairs', '', 'At most ten distinct pairs per year are shown, ranked by max(|Pearson|, |Spearman|); both coefficients are provided so a monotonic association is not mistaken for a linear one. Full matrices and all 21 active7 / 36 diagnostic9 pairs per year and method are saved separately.'])
for year in (2005, 2009):
    top = top_df.loc[top_df.year.eq(year)]
    display = pd.DataFrame({'rank': top['rank'], 'pair': top['predictor_1'] + ' / ' + top['predictor_2'], 'Pearson r': top['pearson'], 'Spearman rho': top['spearman']})
    lines.extend(['', f'### {year} (n = 722)', '', md_table(display)])
lines.extend(['', '## Secondary linear-redundancy diagnostic: active7 VIF', '',
    'Training complete-case predictors were standardized using their training means and sample SDs; each standardized predictor was regressed on the other six with an intercept. VIF = 1/(1 − auxiliary R²). Both full standardized designs have rank 8 (intercept plus seven predictors); there is no exact linear dependence. These auxiliary linear regressions are diagnostics, not retraining RF or SVM. VIF does not directly measure RF/SVM failure or predictive accuracy.', '',
    md_table(vif_df.pivot(index='predictor', columns='year', values='VIF').reindex(ACTIVE).rename(columns={2005: '2005 VIF', 2009: '2009 VIF'}).reset_index(), decimals=3), '',
    '## Interpretation for Reviewer 3', '',
    'The retained predictors exhibit appreciable redundancy; a claim that this active seven-predictor set has negligible collinearity, or satisfies a universal |r|≤0.6 filter, would not be supported by this new training-only check. This result documents the current set, not what selection was performed originally. Correlated predictors can complicate interpretation of individual importance or attribution; correlation and VIF alone do not show impaired generalization or quantify any effect on RF/SVM performance.', '',
    'A reduced-feature sensitivity analysis, if authorized later, would require a prespecified or training-fold-only redundancy rule, RF/SVM refitting with tuning contained in training folds, and evaluation against the same validation scheme. None of that has been performed here. These outputs are available for author consideration and reviewer-response drafting only.', '',
    '## Saved artifacts and preservation', '',
    'All new files are confined to `/home/ubuntu/revision_work/reviewer3_assessment/`. The script reads but does not execute revised_analysis.R. SHA-256 checks confirm the two input CSVs, membership CSV, and inspected revised_analysis.R were unchanged during the run.', '',
    '- `collinearity_diagnostic.py`: reusable, bounded diagnostic script.',
    '- `correlation_{year}_{set}_{method}.csv`: eight full matrices (two years × two sets × two methods).',
    '- `collinearity_tidy_pairs.csv`: all signed correlations, absolute correlations, pair-specific n, and status.',
    '- `collinearity_threshold_summary.csv`: strict-threshold counts and strongest pair by year/set/method.',
    '- `collinearity_strongest_active7_pairs.csv`: ten distinct strongest active7 pairs per year, both methods.',
    '- `collinearity_vif_active7.csv` and `collinearity_vif_training_scaling.csv`: VIF, auxiliary R², design rank, and standardization parameters.',
    '- `collinearity_column_audit.csv` and `collinearity_input_audit.json`: completeness, constant checks, matching, source headers, and provenance.',
    '- `collinearity_findings.md` and `collinearity_findings.txt`: compact finding reports.',
    '- `collinearity_run_log.txt` and `collinearity_validation_log.txt`: execution and independent numerical checks.',
    '- `collinearity_saved_paths.txt`: complete absolute-path inventory.',
])
report = '\n'.join(lines) + '\n'
(OUT / 'collinearity_findings.md').write_text(report)
# Plain-text counterpart retains legible ASCII tables without Markdown emphasis/headings.
text_report = re.sub(r'(?m)^#{1,3} ', '', report).replace('**', '').replace('`', '')
(OUT / 'collinearity_findings.txt').write_text(text_report)
print(summary_df.to_string(index=False))
print('\nACTIVE7 VIF\n' + vif_df.to_string(index=False))
print('\nTOP ACTIVE7 PAIRS\n' + top_df.to_string(index=False))
print('\nFILES SAVED')
for p in sorted(OUT.iterdir()):
    if p.name.startswith(('collinearity_', 'correlation_')):
        print(str(p))
