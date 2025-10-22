/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t prte_mca_oob_tcp_component;

const pmix_mca_base_component_t *prte_oob_base_static_components[] = {
  &prte_mca_oob_tcp_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

