/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_plog_default_component;
extern const pmix_mca_base_component_t pmix_mca_plog_stdfd_component;
extern const pmix_mca_base_component_t pmix_mca_plog_syslog_component;

const pmix_mca_base_component_t *pmix_mca_plog_base_static_components[] = {
  &pmix_mca_plog_default_component, 
  &pmix_mca_plog_stdfd_component, 
  &pmix_mca_plog_syslog_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

