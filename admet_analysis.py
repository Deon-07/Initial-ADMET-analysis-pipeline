#!/usr/bin/env python3
"""
ADMET Analysis Pipeline (final)
Usage: python admet_analysis.py <input_file> [--chunk_size N] [--workers W] ...
"""
from __future__ import annotations

import argparse
import concurrent.futures
import gc
import json
import logging
import math
import sys
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')   # headless backend
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

logger = logging.getLogger("admet_analysis")
logger.setLevel(logging.INFO)
if not logger.handlers:
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(ch)

# Default config
DEFAULT_CONFIG: Dict[str, Any] = {
    'adme_weights': {
        'solubility': 2.0, 'lipophilicity': 1.5, 'absorption': 2.0,
        'bioavailability': 2.0, 'permeability': 1.5, 'cyp_clean': 3.0,
        'bbb': 1.0, 'rule_based': 1.0
    },
    'toxicity_weights': {
        'severe': 5.0, 'moderate': 3.0, 'mild': 1.0, 'pains': 2.0, 'environmental': 1.0
    },
    'thresholds': {
        'logS_pass': -4.0, 'logS_warning': -6.0,
        'logP_optimal_min': 1.0, 'logP_optimal_max': 5.0,
        'logP_warning_min': -1.0, 'logP_warning_max': 7.0
    },
    'hard_failures': {'herg': '+++', 'dili': '+++', 'ames': '+++', 'carcinogenicity': '+++', 'respiratory': '+++'},
    'cyp_weights': {'CYP3A4': 3.0, 'CYP2D6': 2.0, 'CYP2C19': 1.5, 'CYP2C9': 1.0, 'CYP1A2': 1.0},
    'plausibility': {'logP_range': (-5, 10), 'mw_range': (50, 2000)}
}

# Utilities
def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    for k, v in override.items():
        if k in base and isinstance(base[k], dict) and isinstance(v, dict):
            base[k] = deep_merge(base[k], v)
        else:
            base[k] = v
    return base

def safe_float(value: Any, default: float = math.nan, warn: bool = False) -> float:
    try:
        if value is None:
            return default
        if isinstance(value, (float, int, np.integer, np.floating)):
            return float(value)
        s = str(value).strip()
        if s == '':
            return default
        return float(s)
    except Exception:
        return default

def detect_column_patterns(df: pd.DataFrame) -> Dict[str, str]:
    patterns = {
        'solubility': ['logS','LogS','solubility','WaterSol'],
        'lipophilicity': ['logP','LogP','AlogP','MLogP'],
        'absorption': ['HIA','HumanIntestinalAbsorption','Absorption'],
        'permeability': ['Caco2','MDCK','PAMPA','Permeability'],
        'bioavailability': ['F20%','F30%','F(20%)','F(30%)','Bioavailability'],
        'herg': ['hERG','HERG','hERG-inh','Cardiotoxicity'],
        'dili': ['DILI','DrugInducedLiverInjury','Hepatotoxicity'],
        'ames': ['Ames','AMES','Mutagenicity'],
        'pains': ['PAINS','pains','PanAssayInterference']
    }
    detected: Dict[str, str] = {}
    cols = df.columns.tolist()
    for key, candidates in patterns.items():
        for col in cols:
            if any(c.lower() in col.lower() for c in candidates):
                detected[key] = col
                break
    return detected

# Scoring functions adapted for dict inputs
def calculate_adme_score_row(row: Dict[str, Any], df_columns: List[str], config: Dict[str, Any]) -> Tuple[float, Dict[str, Any]]:
    weights = config['adme_weights']; thresholds = config['thresholds']
    score = 0.0; max_possible = 0.0; details: Dict[str, Any] = {}
    # Solubility
    sol_cols = [c for c in df_columns if 'logs' in c.lower() or 'solub' in c.lower()]
    if sol_cols:
        val = safe_float(row.get(sol_cols[0], None))
        if not math.isnan(val):
            if val > thresholds['logS_pass']:
                score += weights['solubility']; details['solubility'] = 'Excellent'
            elif val > thresholds['logS_warning']:
                score += weights['solubility']*0.5; details['solubility'] = 'Acceptable'
            else:
                details['solubility'] = 'Poor'
            max_possible += weights['solubility']
        else:
            details['solubility'] = 'N/A'
    # Lipophilicity
    logp_cols = [c for c in df_columns if 'logp' in c.lower() or 'alogp' in c.lower() or 'mlogp' in c.lower()]
    if logp_cols:
        val = safe_float(row.get(logp_cols[0], None))
        if not math.isnan(val):
            if thresholds['logP_optimal_min'] <= val <= thresholds['logP_optimal_max']:
                score += weights['lipophilicity']; details['lipophilicity'] = 'Optimal'
            elif thresholds['logP_warning_min'] <= val <= thresholds['logP_warning_max']:
                score += weights['lipophilicity'] * 0.5; details['lipophilicity'] = 'Suboptimal'
            else:
                details['lipophilicity'] = 'Poor'
            max_possible += weights['lipophilicity']
        else:
            details['lipophilicity'] = 'N/A'
    # Absorption
    hia_cols = [c for c in df_columns if 'hia' in c.lower() or 'intestinal' in c.lower() or 'absorption' in c.lower()]
    if hia_cols:
        hv = row.get(hia_cols[0], None)
        if hv is not None:
            hia_val = str(hv).strip().lower()
            if 'high' in hia_val or 'good' in hia_val or hia_val in ('yes','true','1'):
                score += weights['absorption']; details['absorption']='Good'
            elif 'medium' in hia_val or 'moderate' in hia_val:
                score += weights['absorption']*0.5; details['absorption']='Moderate'
            else:
                details['absorption']='Poor'
            max_possible += weights['absorption']
    # Bioavailability
    bio_cols = [c for c in df_columns if 'f20' in c.lower() or 'f30' in c.lower() or 'bioavail' in c.lower()]
    if bio_cols:
        bval = row.get(bio_cols[0], None)
        if bval is not None and str(bval).strip().lower() in ('yes','true','1','high'):
            score += weights['bioavailability']; details['bioavailability']='Good'
        else:
            details['bioavailability']='Unknown'
        max_possible += weights['bioavailability']
    # CYP clean profile (approx)
    cyp_weights = config.get('cyp_weights', {})
    cyp_detected=[]; cyp_penalty=False
    for cyp in cyp_weights.keys():
        candidates = [cyp, cyp+'-inh', cyp+'_inhibition', cyp.lower()]
        for cand in df_columns:
            if any(cand.lower()==sc.lower() or sc.lower() in cand.lower() for sc in candidates):
                val=row.get(cand,None)
                if val is not None and str(val).strip().lower() in ('yes','true','1','inhibitor'):
                    cyp_detected.append(cyp); cyp_penalty=True
                break
    if not cyp_penalty:
        score += weights['cyp_clean']; details['cyp_profile']='Clean'
    else:
        details['cyp_profile']=','.join(cyp_detected) if cyp_detected else 'Unknown'
    max_possible += weights['cyp_clean']
    # Rule-based compliance
    rules=['Lipinski','Pfizer','GSK','GoldenTriangle','Veber']; passed=0
    for r in rules:
        if r in df_columns:
            rv=row.get(r,None)
            if rv is not None and str(rv).strip().lower() in ('yes','true','1','pass'):
                passed+=1
    if passed>0:
        rule_weight=(passed/len(rules))*weights['rule_based']; score+=rule_weight; details['rule_based']=f"{passed}/{len(rules)} rules"
    else:
        details['rule_based']='None'
    max_possible+=weights['rule_based']
    normalized=float((score/max_possible*100) if max_possible>0 else 0.0)
    return normalized, details

def calculate_toxicity_penalty_row(row: Dict[str, Any], df_columns: List[str], config: Dict[str, Any]) -> Tuple[float, Dict[str, Any]]:
    weights=config['toxicity_weights']; hard_failures=config['hard_failures']; penalty=0.0; details={}
    severe_list=[]
    for tox,fail_val in hard_failures.items():
        candidates=[c for c in df_columns if tox.lower() in c.lower() or c.lower()==tox.lower()]
        if candidates:
            val=row.get(candidates[0],None)
            if val is not None and str(val).strip()==fail_val:
                penalty+=weights.get('severe',5.0); severe_list.append(tox.upper())
    if severe_list: details['severe_toxicity']=','.join(severe_list)
    moderate_cols=['SkinSen','EyeCorrosion','EyeIrritation','SkinIrritation']; mod_detect=[]
    for col in moderate_cols:
        if col in df_columns:
            v=row.get(col,None)
            if v is not None and str(v).strip().lower() in ('yes','true','1','positive'):
                penalty+=weights.get('moderate',3.0); mod_detect.append(col)
    if mod_detect: details['moderate_toxicity']=','.join(mod_detect)
    nr_cols=['NR-AR','NR-AR-LBD','NR-AhR','NR-Aromatase','NR-ER','NR-ER-LBD','NR-PPAR-gamma','SR-ARE','SR-ATAD5','SR-HSE','SR-MMP','SR-p53']
    mild_detect=[]
    for col in nr_cols:
        if col in df_columns:
            v=row.get(col,None)
            if v is not None and str(v).strip().lower() in ('active','yes','true','1'):
                penalty+=weights.get('mild',1.0); mild_detect.append(col)
    if mild_detect: details['mild_toxicity']=f"{len(mild_detect)} alerts"
    pains_cols=[c for c in df_columns if 'pain' in c.lower()]
    if pains_cols:
        val=row.get(pains_cols[0],None)
        if val is not None and str(val).strip().lower() in ('yes','true','1','alert'):
            penalty+=weights.get('pains',2.0); details['pains']='PAINS'
    return float(penalty), details

def score_record_worker(rec: Dict[str, Any], df_columns: List[str], config: Dict[str, Any]) -> Dict[str, Any]:
    adme_score,adme_detail=calculate_adme_score_row(rec,df_columns,config)
    tox_penalty,tox_detail=calculate_toxicity_penalty_row(rec,df_columns,config)
    final=float(adme_score-tox_penalty)
    return {'ADME_Score':adme_score,'ADME_Details':adme_detail,'Toxicity_Penalty':tox_penalty,'Toxicity_Details':tox_detail,'Final_Score':final}

def batch_score_numeric(df_batch: pd.DataFrame, config: Dict[str, Any]) -> pd.DataFrame:
    thresholds=config['thresholds']; weights=config['adme_weights']; res=pd.DataFrame(index=df_batch.index)
    logp_col=next((c for c in df_batch.columns if 'logp' in c.lower() or 'alogp' in c.lower()), None)
    if logp_col:
        logp_vals=pd.to_numeric(df_batch[logp_col],errors='coerce')
        optimal=logp_vals.between(thresholds['logP_optimal_min'],thresholds['logP_optimal_max'])
        suboptimal=logp_vals.between(thresholds['logP_warning_min'],thresholds['logP_warning_max'])
        lp_pts=np.zeros(len(df_batch),dtype=float)
        lp_pts[optimal.fillna(False)]=weights['lipophilicity']; lp_pts[suboptimal.fillna(False)]=weights['lipophilicity']*0.5
        res['lipophilicity_pts']=lp_pts
    else:
        res['lipophilicity_pts']=0.0
    logs_col=next((c for c in df_batch.columns if 'logs' in c.lower() or 'solub' in c.lower()), None)
    if logs_col:
        logs_vals=pd.to_numeric(df_batch[logs_col],errors='coerce')
        excellent=logs_vals>thresholds['logS_pass']; acceptable=(logs_vals>thresholds['logS_warning']) & (~excellent)
        ls_pts=np.zeros(len(df_batch),dtype=float)
        ls_pts[excellent.fillna(False)]=weights['solubility']; ls_pts[acceptable.fillna(False)]=weights['solubility']*0.5
        res['solubility_pts']=ls_pts
    else:
        res['solubility_pts']=0.0
    res['ADME_Score_partial']=res[['lipophilicity_pts','solubility_pts']].sum(axis=1)
    return res

class ADMETAnalyzer:
    def __init__(self, config: Dict[str, Any] | None = None):
        self.config = deep_merge(DEFAULT_CONFIG.copy(), config or {}) if config else DEFAULT_CONFIG.copy()
        self.df: pd.DataFrame | None = None
        self.results: pd.DataFrame | None = None
        self.detected_cols: Dict[str, str] | None = None

    def load_data(self, input_file: str, sep: str = '\t') -> bool:
        tried=[]
        for delim in [sep,',',';','\t']:
            try:
                df=pd.read_csv(input_file,sep=delim,low_memory=False)
                if df.shape[1]>1:
                    self.df=df; logger.info(f"Loaded {input_file} sep='{delim}' shape={df.shape}"); break
            except Exception as e:
                tried.append((delim,str(e))); continue
        if self.df is None:
            logger.error(f"Failed to load {input_file}. Tried: {tried}"); return False
        self.detected_cols=detect_column_patterns(self.df); logger.info(f"Detected cols: {self.detected_cols}")
        return True

    def check_data_quality(self) -> Dict[str, Any]:
        if self.df is None: raise RuntimeError("No data")
        qc={'rows':len(self.df),'columns':list(self.df.columns),'missing':self.df.isnull().sum().to_dict(),'dtypes':self.df.dtypes.astype(str).to_dict()}
        pl=self.config.get('plausibility',{}); pr={}; logp_col=next((c for c in self.df.columns if 'logp' in c.lower()), None)
        if logp_col:
            logp=pd.to_numeric(self.df[logp_col],errors='coerce'); low,high=pl.get('logP_range',(-5,10)); pr['logP_outliers']=int(((logp<low)|(logp>high)).sum())
        mw_col=next((c for c in self.df.columns if 'molw' in c.lower() or 'mw' in c.lower()), None)
        if mw_col:
            mw=pd.to_numeric(self.df[mw_col],errors='coerce'); low,high=pl.get('mw_range',(50,2000)); pr['molwt_outliers']=int(((mw<low)|(mw>high)).sum())
        qc['plausibility']=pr; logger.info(f"Data quality: {pr}"); return qc

    def _filter_hard_failures(self, df_chunk: pd.DataFrame) -> pd.DataFrame:
        mask=pd.Series(False,index=df_chunk.index)
        for tox,fail_val in self.config['hard_failures'].items():
            candidates=[c for c in df_chunk.columns if tox.lower() in c.lower() or c.lower()==tox.lower()]
            if not candidates: continue
            col=candidates[0]; tox_mask=df_chunk[col].astype(str).str.strip()==fail_val
            if tox_mask.any(): logger.info(f"Filtering removed {tox_mask.sum()} rows where {col}=={fail_val}")
            mask=mask|tox_mask
        return df_chunk.loc[~mask]

    def _expand_details(self, results_df: pd.DataFrame) -> pd.DataFrame:
        def safe_get(d,k):
            try: return d.get(k,'N/A') if isinstance(d,dict) else 'N/A'
            except: return 'N/A'
        results_df['ADME_Solubility']=results_df['ADME_Details'].apply(lambda x: safe_get(x,'solubility'))
        results_df['ADME_Lipophilicity']=results_df['ADME_Details'].apply(lambda x: safe_get(x,'lipophilicity'))
        results_df['ADME_Absorption']=results_df['ADME_Details'].apply(lambda x: safe_get(x,'absorption'))
        results_df['ADME_CYP_Profile']=results_df['ADME_Details'].apply(lambda x: safe_get(x,'cyp_profile'))
        results_df['Toxicity_Severe']=results_df['Toxicity_Details'].apply(lambda x: safe_get(x,'severe_toxicity'))
        results_df['Toxicity_Moderate']=results_df['Toxicity_Details'].apply(lambda x: safe_get(x,'moderate_toxicity'))
        results_df['Toxicity_MildAlerts']=results_df['Toxicity_Details'].apply(lambda x: safe_get(x,'mild_toxicity'))
        results_df['PAINS']=results_df['Toxicity_Details'].apply(lambda x: safe_get(x,'pains'))
        results_df['ADME_Detail_JSON']=results_df['ADME_Details'].apply(lambda x: json.dumps(x) if isinstance(x,dict) else '{}')
        results_df['Tox_Detail_JSON']=results_df['Toxicity_Details'].apply(lambda x: json.dumps(x) if isinstance(x,dict) else '{}')
        for col in ['ADME_Solubility','ADME_Lipophilicity','ADME_Absorption','ADME_CYP_Profile','Toxicity_Severe','PAINS']:
            if col in results_df.columns: results_df[col]=results_df[col].astype('category', copy=False)
        return results_df

    def analyze_compounds(self, top_n:int=20, chunk_size:int=0, parallel_workers:int=0, stream_to_disk:bool=True, output_dir:str='./admet_results') -> pd.DataFrame | None:
        if self.df is None: raise RuntimeError("No data")
        out_path=Path(output_dir); out_path.mkdir(parents=True,exist_ok=True); streamed_csv=out_path/'admet_scored_stream.csv'
        if streamed_csv.exists(): streamed_csv.unlink()
        df_columns=list(self.df.columns); total=len(self.df); logger.info(f"Starting analysis of {total} rows (chunk={chunk_size}, workers={parallel_workers})")
        if chunk_size and chunk_size>0:
            first_write=True
            for start in range(0,total,chunk_size):
                end=min(start+chunk_size,total); chunk=self.df.iloc[start:end]; filtered=self._filter_hard_failures(chunk)
                if filtered.empty: logger.info(f"Chunk {start}:{end} -> 0 survivors"); continue
                records=filtered.to_dict('records')
                try: batch_partial=batch_score_numeric(filtered,self.config)
                except Exception: batch_partial=pd.DataFrame(index=filtered.index)
                if parallel_workers and parallel_workers>1:
                    func=partial(score_record_worker, df_columns=df_columns, config=self.config)
                    with concurrent.futures.ProcessPoolExecutor(max_workers=parallel_workers) as ex:
                        scored_list=list(ex.map(func, records))
                    scored_primitives=pd.DataFrame(scored_list)
                else:
                    scored_primitives=pd.DataFrame([score_record_worker(r,df_columns,self.config) for r in records])
                if 'ADME_Score_partial' in batch_partial.columns:
                    if 'ADME_Score' in scored_primitives.columns:
                        scored_primitives['ADME_Score']=scored_primitives['ADME_Score']+batch_partial['ADME_Score_partial'].values
                        scored_primitives['Final_Score']=scored_primitives['ADME_Score']-scored_primitives['Toxicity_Penalty']
                    else:
                        scored_primitives['ADME_Score']=batch_partial['ADME_Score_partial'].values
                        scored_primitives['Final_Score']=scored_primitives['ADME_Score']-scored_primitives['Toxicity_Penalty']
                merged=pd.concat([filtered.reset_index(drop=True), scored_primitives.reset_index(drop=True)], axis=1)
                merged=self._expand_details(merged)
                if stream_to_disk:
                    merged.to_csv(streamed_csv, mode='a', header=first_write, index=False); first_write=False
                    del chunk, filtered, records, batch_partial, scored_primitives, merged; gc.collect()
                else:
                    if self.results is None: self.results=merged.copy()
                    else: self.results=pd.concat([self.results, merged], ignore_index=True)
                    del chunk, filtered, records, batch_partial, scored_primitives; gc.collect()
            if stream_to_disk: logger.info("Reading back streamed CSV..."); self.results=pd.read_csv(streamed_csv)
        else:
            filtered=self._filter_hard_failures(self.df)
            if filtered.empty: logger.error("No compounds remain after hard failures."); return None
            records=filtered.to_dict('records')
            try: partial_df=batch_score_numeric(filtered, self.config)
            except Exception: partial_df=pd.DataFrame(index=filtered.index)
            if parallel_workers and parallel_workers>1:
                func=partial(score_record_worker, df_columns=df_columns, config=self.config)
                with concurrent.futures.ProcessPoolExecutor(max_workers=parallel_workers) as ex:
                    scored_list=list(ex.map(func, records))
                scored_primitives=pd.DataFrame(scored_list)
            else:
                scored_primitives=pd.DataFrame([score_record_worker(r,df_columns,self.config) for r in records])
            if 'ADME_Score_partial' in partial_df.columns:
                if 'ADME_Score' in scored_primitives.columns:
                    scored_primitives['ADME_Score']=scored_primitives['ADME_Score']+partial_df['ADME_Score_partial'].values
                else:
                    scored_primitives['ADME_Score']=partial_df['ADME_Score_partial'].values
                scored_primitives['Final_Score']=scored_primitives['ADME_Score']-scored_primitives['Toxicity_Penalty']
            self.results=pd.concat([filtered.reset_index(drop=True), scored_primitives.reset_index(drop=True)],axis=1)
            self.results=self._expand_details(self.results)
        if self.results is None or len(self.results)==0: logger.error("No scoring results."); return None
        self.results['Final_Score']=pd.to_numeric(self.results['Final_Score'], errors='coerce').fillna(-999.0)
        self.results['Rank']=self.results['Final_Score'].rank(ascending=False, method='min').astype(int)
        self.results=self.results.sort_values('Rank').reset_index(drop=True)
        logger.info(f"Analysis finished. {len(self.results)} compounds scored. Top score = {self.results['Final_Score'].max():.2f}")
        return self.results.head(top_n)

    def generate_report(self, output_dir: str = './admet_results'):
        if self.results is None: logger.error("No results to report"); return
        out=Path(output_dir); out.mkdir(parents=True, exist_ok=True)
        summary_cols=['Rank','cid','smiles','ADME_Score','Toxicity_Penalty','Final_Score']; available=[c for c in summary_cols if c in self.results.columns]
        self.results[available].to_csv(out/'admet_top_compounds_summary.csv', index=False); self.results.to_csv(out/'admet_comprehensive_results.csv', index=False); self.results.head(20).to_csv(out/'admet_top20_detailed.csv', index=False)
        logger.info(f"Reports saved to {out}")

    def generate_visualizations(self, output_dir: str = './admet_results'):
        if self.results is None: logger.error("No results to visualize"); return
        out=Path(output_dir); out.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(10,7)); sc=plt.scatter(self.results['ADME_Score'], self.results['Toxicity_Penalty'], c=self.results['Final_Score'], cmap='RdYlGn', alpha=0.8, s=40); plt.colorbar(sc,label='Final Score'); plt.xlabel('ADME Score'); plt.ylabel('Toxicity Penalty'); plt.title('ADME vs Toxicity Profile'); plt.grid(alpha=0.25); plt.tight_layout(); plt.savefig(out/'adme_vs_toxicity.png', dpi=300, bbox_inches='tight'); plt.close()
        plt.figure(figsize=(9,6)); plt.hist(self.results['Final_Score'], bins=30, edgecolor='black', alpha=0.7); plt.xlabel('Final Score'); plt.ylabel('Frequency'); plt.title('Final Score Distribution'); plt.grid(axis='y', alpha=0.25); plt.tight_layout(); plt.savefig(out/'final_score_distribution.png', dpi=300, bbox_inches='tight'); plt.close()
        logger.info(f"Visualizations saved to {out}")

    def validate_admet_predictions(self, experimental_df: pd.DataFrame, mapping: Dict[str,str]):
        try:
            from sklearn.metrics import r2_score, mean_absolute_error
            has_sklearn=True
        except Exception:
            has_sklearn=False; logger.warning("scikit-learn not available; fallback to numpy")
        metrics={}
        for pred_col, exp_col in mapping.items():
            if pred_col not in self.results.columns or exp_col not in experimental_df.columns: continue
            y_pred=pd.to_numeric(self.results[pred_col], errors='coerce'); y_true=pd.to_numeric(experimental_df[exp_col], errors='coerce'); mask=y_pred.notna() & y_true.notna(); n=int(mask.sum())
            if n<3: metrics[pred_col]={'n':n,'note':'insufficient data'}; continue
            if has_sklearn:
                r2=float(r2_score(y_true[mask], y_pred[mask])); mae=float(mean_absolute_error(y_true[mask], y_pred[mask]))
            else:
                r=np.corrcoef(y_true[mask], y_pred[mask])[0,1]; r2=float(r**2); mae=float(np.mean(np.abs(y_true[mask]-y_pred[mask])))
            metrics[pred_col]={'n':n,'r2':r2,'mae':mae}
        return metrics

def main():
    parser=argparse.ArgumentParser(description="ADMET analysis pipeline (final)")
    parser.add_argument('input_file', help='Input CSV/TSV file')
    parser.add_argument('-o','--output', default='./admet_results', help='Output dir')
    parser.add_argument('-n','--top_n', type=int, default=20, help='Top N')
    parser.add_argument('--sep', default='\t', help='Input separator')
    parser.add_argument('--config', help='JSON config')
    parser.add_argument('--chunk_size', type=int, default=0, help='Chunk size (0 -> single-shot)')
    parser.add_argument('--workers', type=int, default=0, help='Parallel worker count (0 -> none)')
    parser.add_argument('--no_stream', dest='stream', action='store_false', help='Disable streaming (keep in memory)')
    parser.add_argument('--parquet', action='store_true', help='Write parquet output (pyarrow required)')
    args=parser.parse_args()
    cfg=DEFAULT_CONFIG.copy()
    if args.config:
        try:
            with open(args.config,'r') as fh: userc=json.load(fh); cfg=deep_merge(cfg, userc); logger.info(f"Loaded config {args.config}")
        except Exception as e: logger.warning(f"Failed to load config {args.config}: {e}")
    analyzer=ADMETAnalyzer(config=cfg)
    if not analyzer.load_data(args.input_file, sep=args.sep): logger.error("Load failed"); sys.exit(1)
    analyzer.check_data_quality()
    top_df=analyzer.analyze_compounds(top_n=args.top_n, chunk_size=args.chunk_size, parallel_workers=args.workers, stream_to_disk=args.stream, output_dir=args.output)
    if top_df is None or top_df.empty: logger.error("No results"); sys.exit(1)
    analyzer.generate_report(output_dir=args.output); analyzer.generate_visualizations(output_dir=args.output)
    if args.parquet:
        try: analyzer.results.to_parquet(Path(args.output)/'admet_comprehensive_results.parquet', index=False); logger.info("Wrote parquet")
        except Exception as e: logger.warning(f"Parquet write failed: {e}")
    logger.info("="*60); logger.info("ADMET ANALYSIS SUMMARY"); logger.info(f"Total rows input: {len(analyzer.df)}"); logger.info(f"Rows after scoring: {len(analyzer.results)}"); logger.info(f"Top final score range: {top_df['Final_Score'].min():.2f} - {top_df['Final_Score'].max():.2f}"); logger.info(f"Results saved to: {args.output}"); logger.info("="*60)
    print("\nTop 5 Compounds:")
    for _,r in top_df.head(5).iterrows(): cid=r.get('cid','N/A'); print(f"Rank {int(r['Rank']):2d}: CID={cid:10} | ADME={r['ADME_Score']:6.1f} | Tox={r['Toxicity_Penalty']:6.1f} | Final={r['Final_Score']:7.1f}")

if __name__=='__main__':
    main()
