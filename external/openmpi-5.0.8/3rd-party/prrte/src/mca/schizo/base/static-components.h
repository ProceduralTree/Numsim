/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_schizo_ompi_component;
extern const pmix_mca_base_component_t prte_mca_schizo_prte_component;

const pmix_mca_base_component_t *prte_schizo_base_static_components[] = {
  &prte_mca_schizo_ompi_component, 
  &prte_mca_schizo_prte_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

