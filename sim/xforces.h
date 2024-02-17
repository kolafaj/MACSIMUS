/* some extra forces */
void fixforces(ToIntPtr B, ToIntPtr A);
void centerforces(ToIntPtr B, ToIntPtr A);

void elstforces(ToIntPtr B, ToIntPtr A);

#ifdef SLAB
void slabcutcor(ToIntPtr B, ToIntPtr A);
void wallforces(ToIntPtr B, ToIntPtr A);
#endif

void userforces(ToIntPtr B, ToIntPtr A);

