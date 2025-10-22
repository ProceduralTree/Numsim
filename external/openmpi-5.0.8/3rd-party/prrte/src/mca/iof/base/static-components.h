/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_iof_hnp_component;
extern const pmix_mca_base_component_t prte_mca_iof_prted_component;

const pmix_mca_base_component_t *prte_iof_base_static_components[] = {
  &prte_mca_iof_hnp_component, 
  &prte_mca_iof_prted_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

