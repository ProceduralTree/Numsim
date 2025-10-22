/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const pmix_mca_base_component_t pmix_mca_pinstalldirs_env_component;
extern const pmix_mca_base_component_t pmix_mca_pinstalldirs_config_component;

const pmix_mca_base_component_t *pmix_mca_pinstalldirs_base_static_components[] = {
  &pmix_mca_pinstalldirs_env_component, 
  &pmix_mca_pinstalldirs_config_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

