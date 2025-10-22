/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_ptl_client_component;
extern const pmix_mca_base_component_t pmix_mca_ptl_server_component;
extern const pmix_mca_base_component_t pmix_mca_ptl_tool_component;

const pmix_mca_base_component_t *pmix_mca_ptl_base_static_components[] = {
  &pmix_mca_ptl_client_component, 
  &pmix_mca_ptl_server_component, 
  &pmix_mca_ptl_tool_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

