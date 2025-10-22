/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_psensor_file_component;
extern const pmix_mca_base_component_t pmix_mca_psensor_heartbeat_component;

const pmix_mca_base_component_t *pmix_mca_psensor_base_static_components[] = {
  &pmix_mca_psensor_file_component, 
  &pmix_mca_psensor_heartbeat_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

