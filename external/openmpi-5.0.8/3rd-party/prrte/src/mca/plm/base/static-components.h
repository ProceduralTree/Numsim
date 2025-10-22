/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_plm_slurm_component;
extern const pmix_mca_base_component_t prte_mca_plm_ssh_component;

const pmix_mca_base_component_t *prte_plm_base_static_components[] = {
  &prte_mca_plm_slurm_component, 
  &prte_mca_plm_ssh_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

