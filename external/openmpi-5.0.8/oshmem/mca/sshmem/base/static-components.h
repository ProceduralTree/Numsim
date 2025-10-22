/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const mca_base_component_t mca_sshmem_mmap_component;
extern const mca_base_component_t mca_sshmem_sysv_component;

const mca_base_component_t *mca_sshmem_base_static_components[] = {
  &mca_sshmem_mmap_component, 
  &mca_sshmem_sysv_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

