/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_ess_env_component;
extern const pmix_mca_base_component_t prte_mca_ess_hnp_component;
extern const pmix_mca_base_component_t prte_mca_ess_slurm_component;

const pmix_mca_base_component_t *prte_ess_base_static_components[] = {
  &prte_mca_ess_env_component, 
  &prte_mca_ess_hnp_component, 
  &prte_mca_ess_slurm_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

