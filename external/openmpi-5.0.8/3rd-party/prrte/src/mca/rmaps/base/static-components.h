/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_rmaps_ppr_component;
extern const pmix_mca_base_component_t prte_mca_rmaps_rank_file_component;
extern const pmix_mca_base_component_t prte_mca_rmaps_round_robin_component;
extern const pmix_mca_base_component_t prte_mca_rmaps_seq_component;

const pmix_mca_base_component_t *prte_rmaps_base_static_components[] = {
  &prte_mca_rmaps_ppr_component, 
  &prte_mca_rmaps_rank_file_component, 
  &prte_mca_rmaps_round_robin_component, 
  &prte_mca_rmaps_seq_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

