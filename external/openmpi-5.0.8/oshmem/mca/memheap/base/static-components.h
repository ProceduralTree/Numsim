/*
 * $HEADER$
 */
#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

extern const mca_base_component_t mca_memheap_buddy_component;
extern const mca_base_component_t mca_memheap_ptmalloc_component;

const mca_base_component_t *mca_memheap_base_static_components[] = {
  &mca_memheap_buddy_component, 
  &mca_memheap_ptmalloc_component, 
  NULL
};

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

