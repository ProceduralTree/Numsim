/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_ras_simulator_component;
extern const pmix_mca_base_component_t prte_mca_ras_testrm_component;
extern const pmix_mca_base_component_t prte_mca_ras_pbs_component;
extern const pmix_mca_base_component_t prte_mca_ras_slurm_component;

const pmix_mca_base_component_t *prte_ras_base_static_components[] = {
  &prte_mca_ras_simulator_component, 
  &prte_mca_ras_testrm_component, 
  &prte_mca_ras_pbs_component, 
  &prte_mca_ras_slurm_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

