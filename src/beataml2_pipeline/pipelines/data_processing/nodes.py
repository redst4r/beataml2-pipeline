
import pandas as pd
import numpy as np
import anndata
from pandas.api.types import  is_object_dtype, is_bool_dtype
import mudata

def _subject_id_to_patient(sub: pd.Series):
    # make the patient id a bit more reader-friendly
    return sub.apply(lambda x: f"pt{x}")

def load_mapping(df_sample_mapping: pd.DataFrame):
    # df_sample_mapping = pd.read_excel(BEAT_FOLDER / 'beataml_waves1to4_sample_mapping.xlsx')
    # give the patients some better names than just a number
    df_sample_mapping['dbgap_subject_id'] = _subject_id_to_patient(df_sample_mapping['dbgap_subject_id'])
    return df_sample_mapping

def _consolidate_dna_rna_id_to_labId(rna_series: pd.Series, dna_series:pd.Series, df_mapping:pd.DataFrame):
    """the clinical sheet is missing the labId, but has rnaseq_id and/or dnaseq_id
    reconstruct the labId from those two, using the sample-mapping-datarfame
    """
    rna_to_lab = df_mapping.dropna(subset=['dbgap_rnaseq_sample','labId']).set_index('dbgap_rnaseq_sample').labId.to_dict()
    dna_to_lab = df_mapping.dropna(subset=['dbgap_dnaseq_sample','labId']).set_index('dbgap_dnaseq_sample').labId.to_dict()

    assert len(rna_series) == len(dna_series)

    lab1 = rna_series.apply(lambda x: rna_to_lab[x] if x in rna_to_lab else np.nan)
    lab2 = dna_series.apply(lambda x: dna_to_lab[x] if x in dna_to_lab else np.nan)

    # consolidate
    def _merge_labs(id1, id2):
        if isinstance(id1, str) and isinstance(id2, str):
            assert id1==id2
            return id1
        elif isinstance(id1, str) and not isinstance(id2, str):
            return id1
        elif not isinstance(id1, str) and isinstance(id2, str):
            return id2
        else:
            assert 1==0

    consolidated = [_merge_labs(l1, l2) for l1, l2 in zip(lab1, lab2)]
    return consolidated


def load_clinical(df_clinical: pd.DataFrame, df_mapping: pd.DataFrame):
    # add the labId
    df_clinical['labId'] = _consolidate_dna_rna_id_to_labId(
        df_clinical.dbgap_rnaseq_sample, 
        df_clinical.dbgap_dnaseq_sample, 
        df_mapping)

    # fix the data
    df_clinical['%.Blasts.in.PB'] = df_clinical['%.Blasts.in.PB'].replace({'>90': '90', '>95': '95', 'predominant': np.nan, 'occasional': np.nan, 'rare': np.nan, '75 (by Flow)': np.nan}).astype(float)
    df_clinical['%.Blasts.in.BM']= df_clinical['%.Blasts.in.BM'].replace({'>90': '90', '>95': '95', 'N/A ': np.nan, '<1': 1}).astype(float)
    
    df_clinical['cumulativeTreatmentTypes'] = df_clinical['cumulativeTreatmentTypes'].fillna('unknown')
    df_clinical['had_bm_transplant'] = df_clinical['cumulativeTreatmentTypes'].apply(lambda x: 'y' if 'Bone Marrow Transplant' in x else 'n')
    df_clinical['cumulativeTreatmentStages'] = df_clinical['cumulativeTreatmentStages'].fillna('unknown')
    df_clinical['had_bm_transplant2'] = df_clinical['cumulativeTreatmentStages'].apply(lambda x: 'y' if 'Allogeneic' in x else 'n')

    for c in df_clinical.columns:
        if is_object_dtype(df_clinical[c]):
            df_clinical[c] = df_clinical[c].fillna('unknown')
        # convert bools to int, MuData does not play nice with bools and Nan
        if is_bool_dtype(df_clinical[c]):
            df_clinical[c] = df_clinical[c].astype(int)


    df_clinical = df_clinical.rename({'NPM1': 'NPM1_status', 'RUNX1': 'RUNX1_status', 'ASXL1': 'ASXL1_status','TP53':'TP53_status'}, axis=1)
    return df_clinical


def load_drugs(df_drug, df_family_dict, df_map):
    df_drug_family = df_family_dict["drug_family"]
    df_drug_family_abbr = df_family_dict["synonyms"]
    df_drug_family = df_drug_family.merge(df_drug_family_abbr, on='family')

    # df_drug = pd.read_csv(BEAT_FOLDER / 'beataml_probit_curve_fits_v4_dbgap.txt', sep='\t').query('converged==True')
    df_drug = df_drug.merge(df_drug_family[['inhibitor','family']], on='inhibitor', how='left')
    df_drug['family'] = df_drug['family'].fillna('unknown')
    df_drug['dbgap_subject_id'] = _subject_id_to_patient(df_drug['dbgap_subject_id'])

    df_drug['labId']= _consolidate_dna_rna_id_to_labId(df_drug.dbgap_rnaseq_sample, df_drug.dbgap_dnaseq_sample, df_map)
    df_drug = df_drug.drop(['dbgap_dnaseq_sample','dbgap_rnaseq_sample'], axis=1)


    df_ic50 = pd.crosstab(df_drug.labId, df_drug.inhibitor, values=df_drug.ic50, aggfunc="mean")
    adata_drug_IC50 = anndata.AnnData(df_ic50)

    df_auc = pd.crosstab(df_drug.labId, df_drug.inhibitor, values=df_drug.auc, aggfunc="mean")
    adata_drug_auc = anndata.AnnData(df_auc)

    return adata_drug_IC50, adata_drug_auc


def load_mutations(df_mut, df_mut_suppl, df_mapping, df_clinical_intermediate):

    # rna_to_lab = df_mapping.dropna(subset=['dbgap_rnaseq_sample','labId']).set_index('dbgap_rnaseq_sample').labId.to_dict()
    dna_to_lab = df_mapping.dropna(subset=['dbgap_dnaseq_sample','labId']).set_index('dbgap_dnaseq_sample').labId.to_dict()

    # df_mut_suppl = pd.read_excel(BEAT_FOLDER / 'NIHMS1504008-supplement-Supplementary_Tables_S1-S22.xlsx', sheet_name='Table S7-Variants for Analysis')
    # df_mut = pd.read_csv(BEAT_FOLDER / 'beataml_wes_wv1to4_mutations_dbgap.txt', sep='\t')
    used_hgvsc = set(df_mut_suppl.hgvsc)

    df_mut['labId'] = df_mut['dbgap_sample_id'].apply(lambda x: dna_to_lab[x])
    # filter to used_variants
    df_mut['used_in_manuscript'] = df_mut.hgvsc.isin(used_hgvsc)
    # df_mut= df_mut.query('used_in_manuscript')

    # add patient info
    dna_to_patient = df_mapping.set_index('dbgap_dnaseq_sample').dbgap_subject_id.to_dict()
    df_mut['patient'] = df_mut['dbgap_sample_id'].apply(lambda x: dna_to_patient[x] ) 

    """
    convert to AnnData
    account for the samples that didnt show any mutation and are hence missing from the WES sheet
    """
    df_mut_gene = (pd.crosstab(df_mut.labId, df_mut.symbol) > 0).astype(int)

    # manually add the FLT3 status from the clinical sheet
    flt3_itd_status = df_clinical_intermediate.set_index('labId')['FLT3-ITD'].replace({
        'positive': 1.0,
        'negative': 0.0,
        'unknown': np.nan
    }).to_dict()
    df_mut_gene['FLT3-ITD'] = df_mut_gene.index.map(lambda x: flt3_itd_status[x] if x in flt3_itd_status else np.nan)

    # account for the samples that didnt show any mutation and are hence missing from the WES sheet
    labs_with_dnaseq = list(dna_to_lab.values())
    samples_without_mutation = set(labs_with_dnaseq) - set(df_mut_gene.index)
    for s in samples_without_mutation:
        df_mut_gene.loc[s] = 0

    adata_mut_gene = anndata.AnnData(df_mut_gene)
    adata_mut_gene.X = adata_mut_gene.X.astype(float)
    # adata_mut_gene.X = sparse.csr_matrix(adata_mut_gene.X)
    # adata_mut_gene.X = adata_mut_gene.X.toarray().astype(float)

    return adata_mut_gene

def load_expression(df_expr: pd.DataFrame, df_mapping):
    # df_expr = pd.read_csv(BEAT_FOLDER / 'beataml_waves1to4_norm_exp_dbgap.txt', sep='\t')
    rna_to_lab = df_mapping.dropna(subset=['dbgap_rnaseq_sample','labId']).set_index('dbgap_rnaseq_sample').labId.to_dict()
    # rna_to_patient = df_mapping.set_index('dbgap_rnaseq_sample').dbgap_subject_id.to_dict()
    
    var_names = ['stable_id',	'display_label', 'description','biotype']
    var = df_expr[var_names].set_index('display_label')

    X = df_expr.drop(var_names, axis=1).T
    X.index = X.index.map(lambda x: rna_to_lab[x])

    obs = pd.DataFrame(index=X.index)

    # obs = obs.reset_index().merge(
    #     df_clinical[data_sample_specific], left_on='index', right_on='dbgap_rnaseq_sample', how='left'
    # ).set_index('index').rename({'NPM1': 'NPM1_status', 'RUNX1': 'RUNX1_status', 'ASXL1': 'ASXL1_status','TP53':'TP53_status'}, axis=1)
    # obs['dbgap_subject_id'] = obs.index.map(lambda x: rna_to_patient[x])
    # obs['dbgap_subject_id'] = _subject_id_to_patient(obs['dbgap_subject_id'])
    adata_rna = anndata.AnnData(X, obs=obs, var=var)

    # adata_rna.obs['labId'] = obs.index.map(lambda x: rna_to_lab[x])
    # adata_rna.obs = adata_rna.obs.set_index('labId')

    adata_rna.var.drop('description', axis=1, inplace=True)

    return adata_rna


def load_malignant_celltypes(df: pd.DataFrame, df_mapping: pd.DataFrame):
    # df = catalog.load('malignant_celltypes')
    # df_mapping = catalog.load('mapping_sheet')
    rna_to_lab = df_mapping.dropna(subset=['dbgap_rnaseq_sample','labId']).set_index('dbgap_rnaseq_sample').labId.to_dict()

    df = df.set_index('dbgap_rnaseq_sample')
    df.index = df.index.map(lambda x: rna_to_lab[x])

    return  anndata.AnnData(df)

def to_muon(expr_adata, mut_gene_adata, drug_auc_adata, df_clinical, malignant_ct_adata):
    df_meta = df_clinical.set_index('labId')

    mods = {
        'drugs': drug_auc_adata,
        'rna': expr_adata,
        'dna': mut_gene_adata,
        'malignant_ct': malignant_ct_adata,
    }
    mdata = mudata.MuData(mods)
    mdata.obs = mdata.obs.join(df_meta)
    return mdata

