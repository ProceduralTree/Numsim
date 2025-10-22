/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_errmgr_dvm_component;
extern const pmix_mca_base_component_t prte_mca_errmgr_prted_component;

const pmix_mca_base_component_t *prte_errmgr_base_static_components[] = {
  &prte_mca_errmgr_dvm_component, 
  &prte_mca_errmgr_prted_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

