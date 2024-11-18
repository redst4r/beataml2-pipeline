from kedro.pipeline import Pipeline, node, pipeline

from .nodes import load_expression, load_clinical, load_mapping, load_drugs, load_mutations, to_muon

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=load_mapping,
                inputs="mapping_sheet",
                outputs="mapping_intermediate",
                # name="preprocess_companies_node",
            ),
            node(
                func=load_clinical,
                inputs=["clinical_sheet", "mapping_intermediate"],
                outputs="clinical_intermediate",
                # name="preprocess_companies_node",
            ),
            node(
                func=load_drugs,
                inputs=["drugs", "drug_families","mapping_intermediate"],
                outputs=["drug_ic50_adata", "drug_auc_adata"],
                # name="preprocess_shuttles_node",
            ),
            node(
                func=load_mutations,
                inputs=["mutations", "mutations_suppl", "mapping_intermediate", "clinical_intermediate"],
                outputs="mut_gene_adata",
            ),
            node(
                func=load_expression,
                inputs=["expression", "mapping_intermediate"],
                outputs="expr_adata",
            ),
            node(
                func=to_muon,
                inputs=["expr_adata", "mut_gene_adata", "drug_auc_adata", "clinical_intermediate"],
                outputs="integrated_data",
            ),
        ]
    )
