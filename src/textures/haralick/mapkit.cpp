/* MapKit */

/*
 * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

static char const rcsid[] =
  "@(#) $Jeannot: mapkit.m4c,v 1.23 2004/03/31 19:12:55 js Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include "mapkit.h"
/*
  MapKit
  Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
  @(#) $Jeannot: mapkit.m4,v 1.83 2004/04/10 14:17:03 js Exp $
*/

/* @(#) $Jeannot: mapkit_defs.m4,v 1.20 2004/03/31 19:12:55 js Exp $ */


/* Map structures */

#ifdef MAPKIT_map_ii

/*
  Implementation for map_ii (int -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_ii_keyindex(map_ii *spm, int key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_ii_insertionindex(map_ii *spm, int key);

/* Implementation */

mapkit_error map_ii_init(map_ii *spm)
{
  return map_ii_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_ii_init_hint(map_ii *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0;
  spm->alwaysdefault = 0;

  return map_ii_reallocate(spm, map_ii_meansize(spm, used));
}

mapkit_error map_ii_ensurecapacity(map_ii *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_ii_reallocate(spm, map_ii_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_ii_adjustcapacity(map_ii *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_ii_reallocate(spm, map_ii_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_ii_reallocate(spm, map_ii_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_ii_free(map_ii *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_ii_copy(map_ii *to, map_ii *from)
{
  map_ii_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_ii_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_ii));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_ii_growsize(map_ii *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_ii_shrinksize(map_ii *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_ii_meansize(map_ii *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_ii_reallocate(map_ii *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_ii_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_ii_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    int defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        map_ii_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_ii_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

int map_ii_value_s(map_ii *spm, int key)
{
  mapkit_size_t index;
  
  index = map_ii_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_ii_get_s(map_ii *spm, int key,
  int *value)
{
  mapkit_size_t index;
  
  index = map_ii_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_ii_set_s(map_ii *spm, int key,
  int value)
{
  mapkit_size_t index;

  index = map_ii_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_ii_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_ii_reallocate(spm, map_ii_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

int *map_ii_insertptr_s(map_ii *spm, int key)
{
  mapkit_size_t index;

  index = map_ii_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_ii_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_ii_reallocate(spm, map_ii_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_ii_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

int *map_ii_ptr_s(map_ii *spm, int key)
{
  mapkit_size_t index;

  index = map_ii_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_ii_remove_s(map_ii *spm, int key)
{
  mapkit_size_t index;

  index = map_ii_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_ii_reallocate(spm, map_ii_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_ii_keyindex(map_ii *spm, int key)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_ii_insertionindex(map_ii *spm, int key)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_ii_next(map_ii *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_ii_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_ii_storage *map_ii_nextptr(map_ii *spm, map_ii_storage *pos_contents)
{
  map_ii_storage *end = &(spm->contents[spm->size]);
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_ii_getall(map_ii *spm, map_ii_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_ii_element *pos_array;
  map_ii_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_ii_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_ii_clean(map_ii *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_ii_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_ii_reallocate(spm, map_ii_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_ii_compare(const void *e1, const void *e2)
{
  int key1 = ((map_ii_element *)e1)->key;
  int key2 = ((map_ii_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_ii_getall_sorted(map_ii *spm, map_ii_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_ii_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_ii_compare);

  return MAPKIT_OK;
}

mapkit_error map_ii_setall(map_ii *spm, map_ii_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_ii_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_ii_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_ii_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_ii_removeall(map_ii *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_ii_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_ii_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_ii_printstats(map_ii *spm)
{
  fprintf(stderr, "MAPKIT: map_ii statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_ii */



#ifdef MAPKIT_map_id

/*
  Implementation for map_id (int -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_id_keyindex(map_id *spm, int key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_id_insertionindex(map_id *spm, int key);

/* Implementation */

mapkit_error map_id_init(map_id *spm)
{
  return map_id_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_id_init_hint(map_id *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 0;

  return map_id_reallocate(spm, map_id_meansize(spm, used));
}

mapkit_error map_id_ensurecapacity(map_id *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_id_reallocate(spm, map_id_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_id_adjustcapacity(map_id *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_id_reallocate(spm, map_id_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_id_reallocate(spm, map_id_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_id_free(map_id *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_id_copy(map_id *to, map_id *from)
{
  map_id_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_id_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_id));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_id_growsize(map_id *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_id_shrinksize(map_id *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_id_meansize(map_id *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_id_reallocate(map_id *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_id_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_id_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        map_id_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_id_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double map_id_value_s(map_id *spm, int key)
{
  mapkit_size_t index;
  
  index = map_id_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_id_get_s(map_id *spm, int key,
  double *value)
{
  mapkit_size_t index;
  
  index = map_id_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_id_set_s(map_id *spm, int key,
  double value)
{
  mapkit_size_t index;

  index = map_id_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_id_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_id_reallocate(spm, map_id_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *map_id_insertptr_s(map_id *spm, int key)
{
  mapkit_size_t index;

  index = map_id_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_id_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_id_reallocate(spm, map_id_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_id_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *map_id_ptr_s(map_id *spm, int key)
{
  mapkit_size_t index;

  index = map_id_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_id_remove_s(map_id *spm, int key)
{
  mapkit_size_t index;

  index = map_id_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_id_reallocate(spm, map_id_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_id_keyindex(map_id *spm, int key)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_id_insertionindex(map_id *spm, int key)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_id_next(map_id *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_id_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_id_storage *map_id_nextptr(map_id *spm, map_id_storage *pos_contents)
{
  map_id_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_id_getall(map_id *spm, map_id_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_id_element *pos_array;
  map_id_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_id_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_id_clean(map_id *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_id_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_id_reallocate(spm, map_id_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_id_compare(const void *e1, const void *e2)
{
  int key1 = ((map_id_element *)e1)->key;
  int key2 = ((map_id_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_id_getall_sorted(map_id *spm, map_id_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_id_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_id_compare);

  return MAPKIT_OK;
}

mapkit_error map_id_setall(map_id *spm, map_id_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_id_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_id_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_id_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_id_removeall(map_id *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_id_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_id_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_id_printstats(map_id *spm)
{
  fprintf(stderr, "MAPKIT: map_id statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_id */



#ifdef MAPKIT_map_ivp

/*
  Implementation for map_ivp (int -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_ivp_keyindex(map_ivp *spm, int key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_ivp_insertionindex(map_ivp *spm, int key);

/* Implementation */

mapkit_error map_ivp_init(map_ivp *spm)
{
  return map_ivp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_ivp_init_hint(map_ivp *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = NULL;
  spm->alwaysdefault = 0;

  return map_ivp_reallocate(spm, map_ivp_meansize(spm, used));
}

mapkit_error map_ivp_ensurecapacity(map_ivp *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_ivp_reallocate(spm, map_ivp_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_ivp_adjustcapacity(map_ivp *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_ivp_reallocate(spm, map_ivp_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_ivp_reallocate(spm, map_ivp_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_ivp_free(map_ivp *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_ivp_copy(map_ivp *to, map_ivp *from)
{
  map_ivp_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_ivp_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_ivp));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_ivp_growsize(map_ivp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_ivp_shrinksize(map_ivp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_ivp_meansize(map_ivp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_ivp_reallocate(map_ivp *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_ivp_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_ivp_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    void* defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        map_ivp_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_ivp_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

void* map_ivp_value_s(map_ivp *spm, int key)
{
  mapkit_size_t index;
  
  index = map_ivp_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_ivp_get_s(map_ivp *spm, int key,
  void* *value)
{
  mapkit_size_t index;
  
  index = map_ivp_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_ivp_set_s(map_ivp *spm, int key,
  void* value)
{
  mapkit_size_t index;

  index = map_ivp_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_ivp_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_ivp_reallocate(spm, map_ivp_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

void* *map_ivp_insertptr_s(map_ivp *spm, int key)
{
  mapkit_size_t index;

  index = map_ivp_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_ivp_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_ivp_reallocate(spm, map_ivp_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_ivp_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

void* *map_ivp_ptr_s(map_ivp *spm, int key)
{
  mapkit_size_t index;

  index = map_ivp_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_ivp_remove_s(map_ivp *spm, int key)
{
  mapkit_size_t index;

  index = map_ivp_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_ivp_reallocate(spm, map_ivp_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_ivp_keyindex(map_ivp *spm, int key)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_ivp_insertionindex(map_ivp *spm, int key)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_ivp_next(map_ivp *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_ivp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_ivp_storage *map_ivp_nextptr(map_ivp *spm, map_ivp_storage *pos_contents)
{
  map_ivp_storage *end = &(spm->contents[spm->size]);
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_ivp_getall(map_ivp *spm, map_ivp_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_ivp_element *pos_array;
  map_ivp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_ivp_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_ivp_clean(map_ivp *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_ivp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_ivp_reallocate(spm, map_ivp_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_ivp_compare(const void *e1, const void *e2)
{
  int key1 = ((map_ivp_element *)e1)->key;
  int key2 = ((map_ivp_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_ivp_getall_sorted(map_ivp *spm, map_ivp_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_ivp_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_ivp_compare);

  return MAPKIT_OK;
}

mapkit_error map_ivp_setall(map_ivp *spm, map_ivp_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_ivp_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_ivp_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_ivp_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_ivp_removeall(map_ivp *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_ivp_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_ivp_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_ivp_printstats(map_ivp *spm)
{
  fprintf(stderr, "MAPKIT: map_ivp statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_ivp */




#ifdef MAPKIT_map_h_ii

/*
  Implementation for map_h_ii (int -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_ii_keyindex(map_h_ii *spm, int key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_h_ii_insertionindex(map_h_ii *spm, int key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_ii_init(map_h_ii *spm)
{
  return map_h_ii_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_ii_init_hint(map_h_ii *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0;
  spm->alwaysdefault = 0;

  return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, used));
}

mapkit_error map_h_ii_ensurecapacity(map_h_ii *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_h_ii_adjustcapacity(map_h_ii *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_h_ii_free(map_h_ii *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_ii_copy(map_h_ii *to, map_h_ii *from)
{
  map_h_ii_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_h_ii_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_h_ii));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_h_ii_growsize(map_h_ii *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_ii_shrinksize(map_h_ii *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_h_ii_meansize(map_h_ii *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_ii_reallocate(map_h_ii *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_h_ii_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_h_ii_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    int defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_h_ii_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_h_ii_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

int map_h_ii_value_s(map_h_ii *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_ii_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_h_ii_get_s(map_h_ii *spm, int key,
  int *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_ii_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_h_ii_set_s(map_h_ii *spm, int key,
  int value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ii_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_h_ii_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_h_ii_reallocate(spm, map_h_ii_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

int *map_h_ii_insertptr_s(map_h_ii *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ii_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_h_ii_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_h_ii_reallocate(spm, map_h_ii_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_h_ii_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

int *map_h_ii_ptr_s(map_h_ii *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ii_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_h_ii_remove_s(map_h_ii *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ii_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_ii_reallocate(spm, map_h_ii_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_h_ii_keyindex(map_h_ii *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_h_ii_insertionindex(map_h_ii *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_h_ii_next(map_h_ii *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_h_ii_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_h_ii_storage *map_h_ii_nextptr(map_h_ii *spm, map_h_ii_storage *pos_contents)
{
  map_h_ii_storage *end = &(spm->contents[spm->size]);
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_h_ii_getall(map_h_ii *spm, map_h_ii_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_h_ii_element *pos_array;
  map_h_ii_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_h_ii_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_h_ii_clean(map_h_ii *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_h_ii_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_ii_reallocate(spm, map_h_ii_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_h_ii_compare(const void *e1, const void *e2)
{
  int key1 = ((map_h_ii_element *)e1)->key;
  int key2 = ((map_h_ii_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_ii_getall_sorted(map_h_ii *spm, map_h_ii_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_h_ii_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_h_ii_compare);

  return MAPKIT_OK;
}

mapkit_error map_h_ii_setall(map_h_ii *spm, map_h_ii_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_h_ii_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_ii_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_h_ii_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_h_ii_removeall(map_h_ii *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_ii_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_h_ii_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_h_ii_printstats(map_h_ii *spm)
{
  fprintf(stderr, "MAPKIT: map_h_ii statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_h_ii */



#ifdef MAPKIT_map_h_id

/*
  Implementation for map_h_id (int -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_id_keyindex(map_h_id *spm, int key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_h_id_insertionindex(map_h_id *spm, int key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_id_init(map_h_id *spm)
{
  return map_h_id_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_id_init_hint(map_h_id *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 0;

  return map_h_id_reallocate(spm, map_h_id_meansize(spm, used));
}

mapkit_error map_h_id_ensurecapacity(map_h_id *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_h_id_reallocate(spm, map_h_id_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_h_id_adjustcapacity(map_h_id *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_id_reallocate(spm, map_h_id_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_h_id_reallocate(spm, map_h_id_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_h_id_free(map_h_id *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_id_copy(map_h_id *to, map_h_id *from)
{
  map_h_id_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_h_id_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_h_id));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_h_id_growsize(map_h_id *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_id_shrinksize(map_h_id *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_h_id_meansize(map_h_id *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_id_reallocate(map_h_id *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_h_id_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_h_id_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_h_id_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_h_id_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double map_h_id_value_s(map_h_id *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_id_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_h_id_get_s(map_h_id *spm, int key,
  double *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_id_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_h_id_set_s(map_h_id *spm, int key,
  double value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_id_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_h_id_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_h_id_reallocate(spm, map_h_id_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *map_h_id_insertptr_s(map_h_id *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_id_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_h_id_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_h_id_reallocate(spm, map_h_id_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_h_id_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *map_h_id_ptr_s(map_h_id *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_id_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_h_id_remove_s(map_h_id *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_id_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_id_reallocate(spm, map_h_id_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_h_id_keyindex(map_h_id *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_h_id_insertionindex(map_h_id *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_h_id_next(map_h_id *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_h_id_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_h_id_storage *map_h_id_nextptr(map_h_id *spm, map_h_id_storage *pos_contents)
{
  map_h_id_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_h_id_getall(map_h_id *spm, map_h_id_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_h_id_element *pos_array;
  map_h_id_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_h_id_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_h_id_clean(map_h_id *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_h_id_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_id_reallocate(spm, map_h_id_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_h_id_compare(const void *e1, const void *e2)
{
  int key1 = ((map_h_id_element *)e1)->key;
  int key2 = ((map_h_id_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_id_getall_sorted(map_h_id *spm, map_h_id_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_h_id_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_h_id_compare);

  return MAPKIT_OK;
}

mapkit_error map_h_id_setall(map_h_id *spm, map_h_id_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_h_id_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_id_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_h_id_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_h_id_removeall(map_h_id *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_id_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_h_id_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_h_id_printstats(map_h_id *spm)
{
  fprintf(stderr, "MAPKIT: map_h_id statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_h_id */



#ifdef MAPKIT_map_h_ivp

/*
  Implementation for map_h_ivp (int -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_ivp_keyindex(map_h_ivp *spm, int key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_h_ivp_insertionindex(map_h_ivp *spm, int key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_ivp_init(map_h_ivp *spm)
{
  return map_h_ivp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_ivp_init_hint(map_h_ivp *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = NULL;
  spm->alwaysdefault = 0;

  return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, used));
}

mapkit_error map_h_ivp_ensurecapacity(map_h_ivp *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_h_ivp_adjustcapacity(map_h_ivp *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_h_ivp_free(map_h_ivp *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_ivp_copy(map_h_ivp *to, map_h_ivp *from)
{
  map_h_ivp_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_h_ivp_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_h_ivp));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_h_ivp_growsize(map_h_ivp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_ivp_shrinksize(map_h_ivp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_h_ivp_meansize(map_h_ivp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_ivp_reallocate(map_h_ivp *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_h_ivp_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_h_ivp_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    void* defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_h_ivp_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_h_ivp_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

void* map_h_ivp_value_s(map_h_ivp *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_ivp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_h_ivp_get_s(map_h_ivp *spm, int key,
  void* *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_ivp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_h_ivp_set_s(map_h_ivp *spm, int key,
  void* value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ivp_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_h_ivp_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_h_ivp_reallocate(spm, map_h_ivp_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

void* *map_h_ivp_insertptr_s(map_h_ivp *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ivp_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_h_ivp_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_h_ivp_reallocate(spm, map_h_ivp_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_h_ivp_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

void* *map_h_ivp_ptr_s(map_h_ivp *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ivp_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_h_ivp_remove_s(map_h_ivp *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_ivp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_ivp_reallocate(spm, map_h_ivp_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_h_ivp_keyindex(map_h_ivp *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_h_ivp_insertionindex(map_h_ivp *spm, int key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_h_ivp_next(map_h_ivp *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_h_ivp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_h_ivp_storage *map_h_ivp_nextptr(map_h_ivp *spm, map_h_ivp_storage *pos_contents)
{
  map_h_ivp_storage *end = &(spm->contents[spm->size]);
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_h_ivp_getall(map_h_ivp *spm, map_h_ivp_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_h_ivp_element *pos_array;
  map_h_ivp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_h_ivp_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_h_ivp_clean(map_h_ivp *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_h_ivp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_ivp_reallocate(spm, map_h_ivp_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_h_ivp_compare(const void *e1, const void *e2)
{
  int key1 = ((map_h_ivp_element *)e1)->key;
  int key2 = ((map_h_ivp_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_ivp_getall_sorted(map_h_ivp *spm, map_h_ivp_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_h_ivp_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_h_ivp_compare);

  return MAPKIT_OK;
}

mapkit_error map_h_ivp_setall(map_h_ivp *spm, map_h_ivp_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_h_ivp_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_ivp_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_h_ivp_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_h_ivp_removeall(map_h_ivp *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_ivp_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_h_ivp_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_h_ivp_printstats(map_h_ivp *spm)
{
  fprintf(stderr, "MAPKIT: map_h_ivp statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_h_ivp */




#ifdef MAPKIT_map_vpi

/*
  Implementation for map_vpi (void* -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_vpi_keyindex(map_vpi *spm, void* key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_vpi_insertionindex(map_vpi *spm, void* key);

/* Implementation */

mapkit_error map_vpi_init(map_vpi *spm)
{
  return map_vpi_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_vpi_init_hint(map_vpi *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0;
  spm->alwaysdefault = 0;

  return map_vpi_reallocate(spm, map_vpi_meansize(spm, used));
}

mapkit_error map_vpi_ensurecapacity(map_vpi *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_vpi_reallocate(spm, map_vpi_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_vpi_adjustcapacity(map_vpi *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpi_reallocate(spm, map_vpi_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_vpi_reallocate(spm, map_vpi_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_vpi_free(map_vpi *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_vpi_copy(map_vpi *to, map_vpi *from)
{
  map_vpi_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_vpi_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_vpi));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_vpi_growsize(map_vpi *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_vpi_shrinksize(map_vpi *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_vpi_meansize(map_vpi *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_vpi_reallocate(map_vpi *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_vpi_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_vpi_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    void* key;
    int defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        map_vpi_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_vpi_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

int map_vpi_value_s(map_vpi *spm, void* key)
{
  mapkit_size_t index;
  
  index = map_vpi_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_vpi_get_s(map_vpi *spm, void* key,
  int *value)
{
  mapkit_size_t index;
  
  index = map_vpi_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_vpi_set_s(map_vpi *spm, void* key,
  int value)
{
  mapkit_size_t index;

  index = map_vpi_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_vpi_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_vpi_reallocate(spm, map_vpi_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

int *map_vpi_insertptr_s(map_vpi *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpi_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_vpi_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_vpi_reallocate(spm, map_vpi_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_vpi_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

int *map_vpi_ptr_s(map_vpi *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpi_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_vpi_remove_s(map_vpi *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpi_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpi_reallocate(spm, map_vpi_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_vpi_keyindex(map_vpi *spm, void* key)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_vpi_insertionindex(map_vpi *spm, void* key)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_vpi_next(map_vpi *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_vpi_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_vpi_storage *map_vpi_nextptr(map_vpi *spm, map_vpi_storage *pos_contents)
{
  map_vpi_storage *end = &(spm->contents[spm->size]);
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_vpi_getall(map_vpi *spm, map_vpi_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_vpi_element *pos_array;
  map_vpi_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_vpi_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_vpi_clean(map_vpi *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_vpi_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpi_reallocate(spm, map_vpi_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_vpi_compare(const void *e1, const void *e2)
{
  void* key1 = ((map_vpi_element *)e1)->key;
  void* key2 = ((map_vpi_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_vpi_getall_sorted(map_vpi *spm, map_vpi_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_vpi_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_vpi_compare);

  return MAPKIT_OK;
}

mapkit_error map_vpi_setall(map_vpi *spm, map_vpi_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_vpi_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_vpi_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_vpi_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_vpi_removeall(map_vpi *spm, void* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_vpi_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_vpi_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_vpi_printstats(map_vpi *spm)
{
  fprintf(stderr, "MAPKIT: map_vpi statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_vpi */



#ifdef MAPKIT_map_vpd

/*
  Implementation for map_vpd (void* -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_vpd_keyindex(map_vpd *spm, void* key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_vpd_insertionindex(map_vpd *spm, void* key);

/* Implementation */

mapkit_error map_vpd_init(map_vpd *spm)
{
  return map_vpd_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_vpd_init_hint(map_vpd *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 0;

  return map_vpd_reallocate(spm, map_vpd_meansize(spm, used));
}

mapkit_error map_vpd_ensurecapacity(map_vpd *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_vpd_reallocate(spm, map_vpd_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_vpd_adjustcapacity(map_vpd *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpd_reallocate(spm, map_vpd_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_vpd_reallocate(spm, map_vpd_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_vpd_free(map_vpd *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_vpd_copy(map_vpd *to, map_vpd *from)
{
  map_vpd_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_vpd_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_vpd));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_vpd_growsize(map_vpd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_vpd_shrinksize(map_vpd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_vpd_meansize(map_vpd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_vpd_reallocate(map_vpd *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_vpd_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_vpd_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    void* key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        map_vpd_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_vpd_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double map_vpd_value_s(map_vpd *spm, void* key)
{
  mapkit_size_t index;
  
  index = map_vpd_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_vpd_get_s(map_vpd *spm, void* key,
  double *value)
{
  mapkit_size_t index;
  
  index = map_vpd_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_vpd_set_s(map_vpd *spm, void* key,
  double value)
{
  mapkit_size_t index;

  index = map_vpd_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_vpd_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_vpd_reallocate(spm, map_vpd_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *map_vpd_insertptr_s(map_vpd *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpd_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_vpd_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_vpd_reallocate(spm, map_vpd_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_vpd_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *map_vpd_ptr_s(map_vpd *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpd_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_vpd_remove_s(map_vpd *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpd_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpd_reallocate(spm, map_vpd_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_vpd_keyindex(map_vpd *spm, void* key)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_vpd_insertionindex(map_vpd *spm, void* key)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_vpd_next(map_vpd *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_vpd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_vpd_storage *map_vpd_nextptr(map_vpd *spm, map_vpd_storage *pos_contents)
{
  map_vpd_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_vpd_getall(map_vpd *spm, map_vpd_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_vpd_element *pos_array;
  map_vpd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_vpd_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_vpd_clean(map_vpd *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_vpd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpd_reallocate(spm, map_vpd_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_vpd_compare(const void *e1, const void *e2)
{
  void* key1 = ((map_vpd_element *)e1)->key;
  void* key2 = ((map_vpd_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_vpd_getall_sorted(map_vpd *spm, map_vpd_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_vpd_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_vpd_compare);

  return MAPKIT_OK;
}

mapkit_error map_vpd_setall(map_vpd *spm, map_vpd_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_vpd_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_vpd_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_vpd_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_vpd_removeall(map_vpd *spm, void* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_vpd_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_vpd_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_vpd_printstats(map_vpd *spm)
{
  fprintf(stderr, "MAPKIT: map_vpd statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_vpd */



#ifdef MAPKIT_map_vpvp

/*
  Implementation for map_vpvp (void* -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_vpvp_keyindex(map_vpvp *spm, void* key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_vpvp_insertionindex(map_vpvp *spm, void* key);

/* Implementation */

mapkit_error map_vpvp_init(map_vpvp *spm)
{
  return map_vpvp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_vpvp_init_hint(map_vpvp *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = NULL;
  spm->alwaysdefault = 0;

  return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, used));
}

mapkit_error map_vpvp_ensurecapacity(map_vpvp *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_vpvp_adjustcapacity(map_vpvp *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_vpvp_free(map_vpvp *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_vpvp_copy(map_vpvp *to, map_vpvp *from)
{
  map_vpvp_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_vpvp_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_vpvp));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_vpvp_growsize(map_vpvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_vpvp_shrinksize(map_vpvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_vpvp_meansize(map_vpvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_vpvp_reallocate(map_vpvp *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_vpvp_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_vpvp_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    void* key;
    void* defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        map_vpvp_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_vpvp_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

void* map_vpvp_value_s(map_vpvp *spm, void* key)
{
  mapkit_size_t index;
  
  index = map_vpvp_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_vpvp_get_s(map_vpvp *spm, void* key,
  void* *value)
{
  mapkit_size_t index;
  
  index = map_vpvp_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_vpvp_set_s(map_vpvp *spm, void* key,
  void* value)
{
  mapkit_size_t index;

  index = map_vpvp_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_vpvp_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_vpvp_reallocate(spm, map_vpvp_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

void* *map_vpvp_insertptr_s(map_vpvp *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpvp_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_vpvp_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_vpvp_reallocate(spm, map_vpvp_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_vpvp_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

void* *map_vpvp_ptr_s(map_vpvp *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpvp_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_vpvp_remove_s(map_vpvp *spm, void* key)
{
  mapkit_size_t index;

  index = map_vpvp_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpvp_reallocate(spm, map_vpvp_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_vpvp_keyindex(map_vpvp *spm, void* key)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_vpvp_insertionindex(map_vpvp *spm, void* key)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_vpvp_next(map_vpvp *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_vpvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_vpvp_storage *map_vpvp_nextptr(map_vpvp *spm, map_vpvp_storage *pos_contents)
{
  map_vpvp_storage *end = &(spm->contents[spm->size]);
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_vpvp_getall(map_vpvp *spm, map_vpvp_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_vpvp_element *pos_array;
  map_vpvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_vpvp_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_vpvp_clean(map_vpvp *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_vpvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_vpvp_reallocate(spm, map_vpvp_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_vpvp_compare(const void *e1, const void *e2)
{
  void* key1 = ((map_vpvp_element *)e1)->key;
  void* key2 = ((map_vpvp_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_vpvp_getall_sorted(map_vpvp *spm, map_vpvp_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_vpvp_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_vpvp_compare);

  return MAPKIT_OK;
}

mapkit_error map_vpvp_setall(map_vpvp *spm, map_vpvp_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_vpvp_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_vpvp_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_vpvp_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_vpvp_removeall(map_vpvp *spm, void* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_vpvp_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_vpvp_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_vpvp_printstats(map_vpvp *spm)
{
  fprintf(stderr, "MAPKIT: map_vpvp statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_vpvp */




#ifdef MAPKIT_map_h_vpi

/*
  Implementation for map_h_vpi (void* -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_vpi_keyindex(map_h_vpi *spm, void* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_h_vpi_insertionindex(map_h_vpi *spm, void* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_vpi_init(map_h_vpi *spm)
{
  return map_h_vpi_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_vpi_init_hint(map_h_vpi *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0;
  spm->alwaysdefault = 0;

  return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, used));
}

mapkit_error map_h_vpi_ensurecapacity(map_h_vpi *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_h_vpi_adjustcapacity(map_h_vpi *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_h_vpi_free(map_h_vpi *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_vpi_copy(map_h_vpi *to, map_h_vpi *from)
{
  map_h_vpi_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_h_vpi_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_h_vpi));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_h_vpi_growsize(map_h_vpi *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_vpi_shrinksize(map_h_vpi *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_h_vpi_meansize(map_h_vpi *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_vpi_reallocate(map_h_vpi *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_h_vpi_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_h_vpi_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    void* key;
    int defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_h_vpi_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_h_vpi_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

int map_h_vpi_value_s(map_h_vpi *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_vpi_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_h_vpi_get_s(map_h_vpi *spm, void* key,
  int *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_vpi_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_h_vpi_set_s(map_h_vpi *spm, void* key,
  int value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpi_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_h_vpi_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_h_vpi_reallocate(spm, map_h_vpi_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

int *map_h_vpi_insertptr_s(map_h_vpi *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpi_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_h_vpi_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_h_vpi_reallocate(spm, map_h_vpi_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_h_vpi_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

int *map_h_vpi_ptr_s(map_h_vpi *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpi_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_h_vpi_remove_s(map_h_vpi *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpi_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpi_reallocate(spm, map_h_vpi_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_h_vpi_keyindex(map_h_vpi *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_h_vpi_insertionindex(map_h_vpi *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_h_vpi_next(map_h_vpi *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_h_vpi_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_h_vpi_storage *map_h_vpi_nextptr(map_h_vpi *spm, map_h_vpi_storage *pos_contents)
{
  map_h_vpi_storage *end = &(spm->contents[spm->size]);
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_h_vpi_getall(map_h_vpi *spm, map_h_vpi_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_h_vpi_element *pos_array;
  map_h_vpi_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_h_vpi_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_h_vpi_clean(map_h_vpi *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_h_vpi_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpi_reallocate(spm, map_h_vpi_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_h_vpi_compare(const void *e1, const void *e2)
{
  void* key1 = ((map_h_vpi_element *)e1)->key;
  void* key2 = ((map_h_vpi_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_vpi_getall_sorted(map_h_vpi *spm, map_h_vpi_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_h_vpi_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_h_vpi_compare);

  return MAPKIT_OK;
}

mapkit_error map_h_vpi_setall(map_h_vpi *spm, map_h_vpi_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_h_vpi_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_vpi_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_h_vpi_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_h_vpi_removeall(map_h_vpi *spm, void* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_vpi_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_h_vpi_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_h_vpi_printstats(map_h_vpi *spm)
{
  fprintf(stderr, "MAPKIT: map_h_vpi statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_h_vpi */



#ifdef MAPKIT_map_h_vpd

/*
  Implementation for map_h_vpd (void* -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_vpd_keyindex(map_h_vpd *spm, void* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_h_vpd_insertionindex(map_h_vpd *spm, void* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_vpd_init(map_h_vpd *spm)
{
  return map_h_vpd_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_vpd_init_hint(map_h_vpd *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 0;

  return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, used));
}

mapkit_error map_h_vpd_ensurecapacity(map_h_vpd *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_h_vpd_adjustcapacity(map_h_vpd *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_h_vpd_free(map_h_vpd *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_vpd_copy(map_h_vpd *to, map_h_vpd *from)
{
  map_h_vpd_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_h_vpd_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_h_vpd));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_h_vpd_growsize(map_h_vpd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_vpd_shrinksize(map_h_vpd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_h_vpd_meansize(map_h_vpd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_vpd_reallocate(map_h_vpd *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_h_vpd_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_h_vpd_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    void* key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_h_vpd_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_h_vpd_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double map_h_vpd_value_s(map_h_vpd *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_vpd_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_h_vpd_get_s(map_h_vpd *spm, void* key,
  double *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_vpd_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_h_vpd_set_s(map_h_vpd *spm, void* key,
  double value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpd_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_h_vpd_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_h_vpd_reallocate(spm, map_h_vpd_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *map_h_vpd_insertptr_s(map_h_vpd *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpd_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_h_vpd_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_h_vpd_reallocate(spm, map_h_vpd_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_h_vpd_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *map_h_vpd_ptr_s(map_h_vpd *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpd_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_h_vpd_remove_s(map_h_vpd *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpd_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpd_reallocate(spm, map_h_vpd_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_h_vpd_keyindex(map_h_vpd *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_h_vpd_insertionindex(map_h_vpd *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_h_vpd_next(map_h_vpd *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_h_vpd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_h_vpd_storage *map_h_vpd_nextptr(map_h_vpd *spm, map_h_vpd_storage *pos_contents)
{
  map_h_vpd_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_h_vpd_getall(map_h_vpd *spm, map_h_vpd_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_h_vpd_element *pos_array;
  map_h_vpd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_h_vpd_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_h_vpd_clean(map_h_vpd *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_h_vpd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpd_reallocate(spm, map_h_vpd_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_h_vpd_compare(const void *e1, const void *e2)
{
  void* key1 = ((map_h_vpd_element *)e1)->key;
  void* key2 = ((map_h_vpd_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_vpd_getall_sorted(map_h_vpd *spm, map_h_vpd_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_h_vpd_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_h_vpd_compare);

  return MAPKIT_OK;
}

mapkit_error map_h_vpd_setall(map_h_vpd *spm, map_h_vpd_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_h_vpd_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_vpd_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_h_vpd_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_h_vpd_removeall(map_h_vpd *spm, void* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_vpd_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_h_vpd_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_h_vpd_printstats(map_h_vpd *spm)
{
  fprintf(stderr, "MAPKIT: map_h_vpd statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_h_vpd */



#ifdef MAPKIT_map_h_vpvp

/*
  Implementation for map_h_vpvp (void* -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_h_vpvp_keyindex(map_h_vpvp *spm, void* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_h_vpvp_insertionindex(map_h_vpvp *spm, void* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_h_vpvp_init(map_h_vpvp *spm)
{
  return map_h_vpvp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_h_vpvp_init_hint(map_h_vpvp *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = NULL;
  spm->alwaysdefault = 0;

  return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, used));
}

mapkit_error map_h_vpvp_ensurecapacity(map_h_vpvp *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_h_vpvp_adjustcapacity(map_h_vpvp *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_h_vpvp_free(map_h_vpvp *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_h_vpvp_copy(map_h_vpvp *to, map_h_vpvp *from)
{
  map_h_vpvp_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_h_vpvp_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_h_vpvp));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_h_vpvp_growsize(map_h_vpvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_h_vpvp_shrinksize(map_h_vpvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_h_vpvp_meansize(map_h_vpvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_h_vpvp_reallocate(map_h_vpvp *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_h_vpvp_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_h_vpvp_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    void* key;
    void* defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_h_vpvp_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_hash((mapkit_hash_t) key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_h_vpvp_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

void* map_h_vpvp_value_s(map_h_vpvp *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_vpvp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_h_vpvp_get_s(map_h_vpvp *spm, void* key,
  void* *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_h_vpvp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_h_vpvp_set_s(map_h_vpvp *spm, void* key,
  void* value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpvp_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_h_vpvp_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_h_vpvp_reallocate(spm, map_h_vpvp_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

void* *map_h_vpvp_insertptr_s(map_h_vpvp *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpvp_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_h_vpvp_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_h_vpvp_reallocate(spm, map_h_vpvp_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_h_vpvp_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

void* *map_h_vpvp_ptr_s(map_h_vpvp *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpvp_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_h_vpvp_remove_s(map_h_vpvp *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_h_vpvp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpvp_reallocate(spm, map_h_vpvp_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_h_vpvp_keyindex(map_h_vpvp *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! ((spm->contents[index].key) == (key)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_h_vpvp_insertionindex(map_h_vpvp *spm, void* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && ((spm->contents[index].key) == (key))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! ((spm->contents[index].key) == (key))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! ((spm->contents[index].key) == (key)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_h_vpvp_next(map_h_vpvp *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_h_vpvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_h_vpvp_storage *map_h_vpvp_nextptr(map_h_vpvp *spm, map_h_vpvp_storage *pos_contents)
{
  map_h_vpvp_storage *end = &(spm->contents[spm->size]);
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_h_vpvp_getall(map_h_vpvp *spm, map_h_vpvp_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_h_vpvp_element *pos_array;
  map_h_vpvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_h_vpvp_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_h_vpvp_clean(map_h_vpvp *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_h_vpvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_h_vpvp_reallocate(spm, map_h_vpvp_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_h_vpvp_compare(const void *e1, const void *e2)
{
  void* key1 = ((map_h_vpvp_element *)e1)->key;
  void* key2 = ((map_h_vpvp_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error map_h_vpvp_getall_sorted(map_h_vpvp *spm, map_h_vpvp_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_h_vpvp_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_h_vpvp_compare);

  return MAPKIT_OK;
}

mapkit_error map_h_vpvp_setall(map_h_vpvp *spm, map_h_vpvp_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_h_vpvp_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_vpvp_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_h_vpvp_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_h_vpvp_removeall(map_h_vpvp *spm, void* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_h_vpvp_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_h_vpvp_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_h_vpvp_printstats(map_h_vpvp *spm)
{
  fprintf(stderr, "MAPKIT: map_h_vpvp statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_h_vpvp */




#ifdef MAPKIT_map_stri

/*
  Implementation for map_stri (char* -> int)
  Default value : 0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_stri_keyindex(map_stri *spm, char* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_stri_insertionindex(map_stri *spm, char* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_stri_init(map_stri *spm)
{
  return map_stri_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_stri_init_hint(map_stri *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0;
  spm->alwaysdefault = 0;

  return map_stri_reallocate(spm, map_stri_meansize(spm, used));
}

mapkit_error map_stri_ensurecapacity(map_stri *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_stri_reallocate(spm, map_stri_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_stri_adjustcapacity(map_stri *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_stri_reallocate(spm, map_stri_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_stri_reallocate(spm, map_stri_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_stri_free(map_stri *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_stri_copy(map_stri *to, map_stri *from)
{
  map_stri_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_stri_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_stri));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_stri_growsize(map_stri *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_stri_shrinksize(map_stri *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_stri_meansize(map_stri *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_stri_reallocate(map_stri *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_stri_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_stri_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    char* key;
    int defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_stri_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_strhash(key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_stri_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

int map_stri_value_s(map_stri *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_stri_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_stri_get_s(map_stri *spm, char* key,
  int *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_stri_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_stri_set_s(map_stri *spm, char* key,
  int value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_stri_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_stri_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_stri_reallocate(spm, map_stri_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

int *map_stri_insertptr_s(map_stri *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_stri_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_stri_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_stri_reallocate(spm, map_stri_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_stri_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

int *map_stri_ptr_s(map_stri *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_stri_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_stri_remove_s(map_stri *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_stri_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_stri_reallocate(spm, map_stri_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_stri_keyindex(map_stri *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! (strcmp(spm->contents[index].key,key) == 0))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_stri_insertionindex(map_stri *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && (strcmp(spm->contents[index].key,key) == 0)) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! (strcmp(spm->contents[index].key,key) == 0)))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! (strcmp(spm->contents[index].key,key) == 0))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_stri_next(map_stri *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_stri_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_stri_storage *map_stri_nextptr(map_stri *spm, map_stri_storage *pos_contents)
{
  map_stri_storage *end = &(spm->contents[spm->size]);
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_stri_getall(map_stri *spm, map_stri_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_stri_element *pos_array;
  map_stri_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_stri_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_stri_clean(map_stri *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_stri_storage *pos_contents;
  int defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_stri_reallocate(spm, map_stri_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_stri_compare(const void *e1, const void *e2)
{
  char* key1 = ((map_stri_element *)e1)->key;
  char* key2 = ((map_stri_element *)e2)->key;

  return strcmp(key1,key2);
}

mapkit_error map_stri_getall_sorted(map_stri *spm, map_stri_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_stri_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_stri_compare);

  return MAPKIT_OK;
}

mapkit_error map_stri_setall(map_stri *spm, map_stri_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_stri_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_stri_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_stri_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_stri_removeall(map_stri *spm, char* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_stri_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_stri_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_stri_printstats(map_stri *spm)
{
  fprintf(stderr, "MAPKIT: map_stri statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_stri */



#ifdef MAPKIT_map_strd

/*
  Implementation for map_strd (char* -> double)
  Default value : 0.0
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_strd_keyindex(map_strd *spm, char* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_strd_insertionindex(map_strd *spm, char* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_strd_init(map_strd *spm)
{
  return map_strd_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_strd_init_hint(map_strd *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 0;

  return map_strd_reallocate(spm, map_strd_meansize(spm, used));
}

mapkit_error map_strd_ensurecapacity(map_strd *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_strd_reallocate(spm, map_strd_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_strd_adjustcapacity(map_strd *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strd_reallocate(spm, map_strd_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_strd_reallocate(spm, map_strd_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_strd_free(map_strd *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_strd_copy(map_strd *to, map_strd *from)
{
  map_strd_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_strd_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_strd));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_strd_growsize(map_strd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_strd_shrinksize(map_strd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_strd_meansize(map_strd *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_strd_reallocate(map_strd *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_strd_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_strd_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    char* key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_strd_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_strhash(key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_strd_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double map_strd_value_s(map_strd *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_strd_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_strd_get_s(map_strd *spm, char* key,
  double *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_strd_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_strd_set_s(map_strd *spm, char* key,
  double value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strd_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_strd_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_strd_reallocate(spm, map_strd_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *map_strd_insertptr_s(map_strd *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strd_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_strd_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_strd_reallocate(spm, map_strd_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_strd_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *map_strd_ptr_s(map_strd *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strd_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_strd_remove_s(map_strd *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strd_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strd_reallocate(spm, map_strd_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_strd_keyindex(map_strd *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! (strcmp(spm->contents[index].key,key) == 0))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_strd_insertionindex(map_strd *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && (strcmp(spm->contents[index].key,key) == 0)) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! (strcmp(spm->contents[index].key,key) == 0)))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! (strcmp(spm->contents[index].key,key) == 0))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_strd_next(map_strd *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_strd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_strd_storage *map_strd_nextptr(map_strd *spm, map_strd_storage *pos_contents)
{
  map_strd_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_strd_getall(map_strd *spm, map_strd_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_strd_element *pos_array;
  map_strd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_strd_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_strd_clean(map_strd *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_strd_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strd_reallocate(spm, map_strd_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_strd_compare(const void *e1, const void *e2)
{
  char* key1 = ((map_strd_element *)e1)->key;
  char* key2 = ((map_strd_element *)e2)->key;

  return strcmp(key1,key2);
}

mapkit_error map_strd_getall_sorted(map_strd *spm, map_strd_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_strd_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_strd_compare);

  return MAPKIT_OK;
}

mapkit_error map_strd_setall(map_strd *spm, map_strd_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_strd_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_strd_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_strd_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_strd_removeall(map_strd *spm, char* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_strd_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_strd_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_strd_printstats(map_strd *spm)
{
  fprintf(stderr, "MAPKIT: map_strd statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_strd */



#ifdef MAPKIT_map_strvp

/*
  Implementation for map_strvp (char* -> void*)
  Default value : NULL
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_strvp_keyindex(map_strvp *spm, char* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_strvp_insertionindex(map_strvp *spm, char* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_strvp_init(map_strvp *spm)
{
  return map_strvp_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_strvp_init_hint(map_strvp *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = NULL;
  spm->alwaysdefault = 0;

  return map_strvp_reallocate(spm, map_strvp_meansize(spm, used));
}

mapkit_error map_strvp_ensurecapacity(map_strvp *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_strvp_reallocate(spm, map_strvp_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_strvp_adjustcapacity(map_strvp *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strvp_reallocate(spm, map_strvp_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_strvp_reallocate(spm, map_strvp_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_strvp_free(map_strvp *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_strvp_copy(map_strvp *to, map_strvp *from)
{
  map_strvp_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_strvp_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_strvp));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_strvp_growsize(map_strvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_strvp_shrinksize(map_strvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_strvp_meansize(map_strvp *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_strvp_reallocate(map_strvp *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_strvp_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_strvp_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    char* key;
    void* defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_strvp_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_strhash(key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_strvp_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

void* map_strvp_value_s(map_strvp *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_strvp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_strvp_get_s(map_strvp *spm, char* key,
  void* *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_strvp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_strvp_set_s(map_strvp *spm, char* key,
  void* value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strvp_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_strvp_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_strvp_reallocate(spm, map_strvp_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

void* *map_strvp_insertptr_s(map_strvp *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strvp_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_strvp_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_strvp_reallocate(spm, map_strvp_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_strvp_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

void* *map_strvp_ptr_s(map_strvp *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strvp_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_strvp_remove_s(map_strvp *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strvp_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strvp_reallocate(spm, map_strvp_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_strvp_keyindex(map_strvp *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! (strcmp(spm->contents[index].key,key) == 0))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_strvp_insertionindex(map_strvp *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && (strcmp(spm->contents[index].key,key) == 0)) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! (strcmp(spm->contents[index].key,key) == 0)))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! (strcmp(spm->contents[index].key,key) == 0))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_strvp_next(map_strvp *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_strvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

map_strvp_storage *map_strvp_nextptr(map_strvp *spm, map_strvp_storage *pos_contents)
{
  map_strvp_storage *end = &(spm->contents[spm->size]);
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error map_strvp_getall(map_strvp *spm, map_strvp_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_strvp_element *pos_array;
  map_strvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_strvp_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_strvp_clean(map_strvp *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_strvp_storage *pos_contents;
  void* defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strvp_reallocate(spm, map_strvp_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_strvp_compare(const void *e1, const void *e2)
{
  char* key1 = ((map_strvp_element *)e1)->key;
  char* key2 = ((map_strvp_element *)e2)->key;

  return strcmp(key1,key2);
}

mapkit_error map_strvp_getall_sorted(map_strvp *spm, map_strvp_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_strvp_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_strvp_compare);

  return MAPKIT_OK;
}

mapkit_error map_strvp_setall(map_strvp *spm, map_strvp_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_strvp_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_strvp_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_strvp_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_strvp_removeall(map_strvp *spm, char* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_strvp_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_strvp_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_strvp_printstats(map_strvp *spm)
{
  fprintf(stderr, "MAPKIT: map_strvp statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_strvp */



#ifdef MAPKIT_map_strstr

/*
  Implementation for map_strstr (char* -> char*)
  Default value : ""
  Uses a state field.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t map_strstr_keyindex(map_strstr *spm, char* key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t map_strstr_insertionindex(map_strstr *spm, char* key, mapkit_hash_t hash);

/* Implementation */

mapkit_error map_strstr_init(map_strstr *spm)
{
  return map_strstr_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error map_strstr_init_hint(map_strstr *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = (char *)"";
  spm->alwaysdefault = 0;

  return map_strstr_reallocate(spm, map_strstr_meansize(spm, used));
}

mapkit_error map_strstr_ensurecapacity(map_strstr *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return map_strstr_reallocate(spm, map_strstr_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error map_strstr_adjustcapacity(map_strstr *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strstr_reallocate(spm, map_strstr_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return map_strstr_reallocate(spm, map_strstr_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void map_strstr_free(map_strstr *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error map_strstr_copy(map_strstr *to, map_strstr *from)
{
  map_strstr_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (map_strstr_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(map_strstr));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t map_strstr_growsize(map_strstr *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t map_strstr_shrinksize(map_strstr *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t map_strstr_meansize(map_strstr *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error map_strstr_reallocate(map_strstr *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  map_strstr_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (map_strstr_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    char* key;
    char* defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        map_strstr_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = mapkit_strhash(key)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = map_strstr_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! (strcmp(oldcontents[index].value,defaultvalue) == 0)))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

char* map_strstr_value_s(map_strstr *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_strstr_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error map_strstr_get_s(map_strstr *spm, char* key,
  char* *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = map_strstr_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error map_strstr_set_s(map_strstr *spm, char* key,
  char* value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strstr_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    map_strstr_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return map_strstr_reallocate(spm, map_strstr_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

char* *map_strstr_insertptr_s(map_strstr *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strstr_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    map_strstr_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = map_strstr_reallocate(spm, map_strstr_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = map_strstr_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

char* *map_strstr_ptr_s(map_strstr *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strstr_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! (strcmp(spm->contents[index].value,spm->defaultvalue) == 0))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error map_strstr_remove_s(map_strstr *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = map_strstr_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strstr_reallocate(spm, map_strstr_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t map_strstr_keyindex(map_strstr *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! (strcmp(spm->contents[index].key,key) == 0))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t map_strstr_insertionindex(map_strstr *spm, char* key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && (strcmp(spm->contents[index].key,key) == 0)) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! (strcmp(spm->contents[index].key,key) == 0)))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! (strcmp(spm->contents[index].key,key) == 0))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t map_strstr_next(map_strstr *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  map_strstr_storage *pos_contents;
  char* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! (strcmp(pos_contents->value,defaultvalue) == 0))))
      return index;

  return -1;
}

map_strstr_storage *map_strstr_nextptr(map_strstr *spm, map_strstr_storage *pos_contents)
{
  map_strstr_storage *end = &(spm->contents[spm->size]);
  char* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! (strcmp(pos_contents->value,defaultvalue) == 0))))
      return pos_contents;

  return NULL;
}

mapkit_error map_strstr_getall(map_strstr *spm, map_strstr_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  map_strstr_element *pos_array;
  map_strstr_storage *pos_contents;
  char* defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (map_strstr_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! (strcmp(pos_contents->value,defaultvalue) == 0))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error map_strstr_clean(map_strstr *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  map_strstr_storage *pos_contents;
  char* defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (strcmp(pos_contents->value,defaultvalue) == 0))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return map_strstr_reallocate(spm, map_strstr_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int map_strstr_compare(const void *e1, const void *e2)
{
  char* key1 = ((map_strstr_element *)e1)->key;
  char* key2 = ((map_strstr_element *)e2)->key;

  return strcmp(key1,key2);
}

mapkit_error map_strstr_getall_sorted(map_strstr *spm, map_strstr_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = map_strstr_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), map_strstr_compare);

  return MAPKIT_OK;
}

mapkit_error map_strstr_setall(map_strstr *spm, map_strstr_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = map_strstr_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_strstr_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    map_strstr_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error map_strstr_removeall(map_strstr *spm, char* *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = map_strstr_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  map_strstr_adjustcapacity(spm);

  return MAPKIT_OK;
}

void map_strstr_printstats(map_strstr *spm)
{
  fprintf(stderr, "MAPKIT: map_strstr statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_map_strstr */



/* Sparse structures */

#ifdef MAPKIT_spvector

/*
  Implementation for spvector (int -> double)
  Default value : 0.0
  Unassigned values always default to the default value.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t spvector_keyindex(spvector *spm, int key);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t spvector_insertionindex(spvector *spm, int key);

/* Implementation */

mapkit_error spvector_init(spvector *spm)
{
  return spvector_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error spvector_init_hint(spvector *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 1;

  return spvector_reallocate(spm, spvector_meansize(spm, used));
}

mapkit_error spvector_ensurecapacity(spvector *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return spvector_reallocate(spm, spvector_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error spvector_adjustcapacity(spvector *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return spvector_reallocate(spm, spvector_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return spvector_reallocate(spm, spvector_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void spvector_free(spvector *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error spvector_copy(spvector *to, spvector *from)
{
  spvector_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (spvector_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(spvector));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t spvector_growsize(spvector *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t spvector_shrinksize(spvector *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t spvector_meansize(spvector *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error spvector_reallocate(spvector *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  spvector_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = (spvector_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].key = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    int key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if ((key = oldcontents[index].key) >= MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        spvector_storage *contents;

        /* Fast path */
        ins_index = ((mapkit_hash_t)key) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->key != MAPKIT_FREESLOT)
        {
          ins_index = spvector_insertionindex(spm, key);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double spvector_value_s(spvector *spm, int key)
{
  mapkit_size_t index;
  
  index = spvector_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error spvector_get_s(spvector *spm, int key,
  double *value)
{
  mapkit_size_t index;
  
  index = spvector_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error spvector_set_s(spvector *spm, int key,
  double value)
{
  mapkit_size_t index;

  index = spvector_insertionindex(spm, key);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    spvector_storage *element = &(spm->contents[index]);
    int free = element->key == MAPKIT_FREESLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return spvector_reallocate(spm, spvector_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *spvector_insertptr_s(spvector *spm, int key)
{
  mapkit_size_t index;

  index = spvector_insertionindex(spm, key);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    spvector_storage *element = &(spm->contents[index]);
    if (element->key == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = spvector_reallocate(spm, spvector_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = spvector_insertionindex(spm, key);
        /* FREESLOT */
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *spvector_ptr_s(spvector *spm, int key)
{
  mapkit_size_t index;

  index = spvector_keyindex(spm, key);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error spvector_remove_s(spvector *spm, int key)
{
  mapkit_size_t index;

  index = spvector_keyindex(spm, key);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].key = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return spvector_reallocate(spm, spvector_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t spvector_keyindex(spvector *spm, int key)
{
  mapkit_size_t index, decrement;

  int ckey;

#ifdef MAPKIT_DEBUG
  if (key < MAPKIT_FULLSLOT)
  {
    /*   Not user-called : should never happen */
    MAPKIT_ERROR_NORET(MAPKIT_EBADKEY);
    return MAPKIT_KEYNOTFOUND;
  }
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;
  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((ckey = spm->contents[index].key) != MAPKIT_FREESLOT
    && (ckey == MAPKIT_DELETEDSLOT || ckey != key))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (ckey == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t spvector_insertionindex(spvector *spm, int key)
{
  mapkit_size_t index, decrement;
  int ckey;

#ifdef MAPKIT_DEBUG
  if (key < MAPKIT_FULLSLOT)
    /*   Not user-called : should never happen */
    MAPKIT_FATAL_ERROR(MAPKIT_EBADKEY);
#endif

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = ((mapkit_hash_t)key) % spm->size;

  /* Fast path (largely superfluous) */
  if ((ckey = spm->contents[index].key) == MAPKIT_FREESLOT) return index;
  if (ckey == key) return -index-1;

  decrement = (((mapkit_hash_t)key) % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((ckey >= 0) && (ckey != key))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    ckey = spm->contents[index].key;
  }

  if (ckey == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((ckey != MAPKIT_FREESLOT)
      && (ckey == MAPKIT_DELETEDSLOT || ckey != key))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      ckey = spm->contents[index].key;
    }
    if (ckey == MAPKIT_FREESLOT) return index2;
  }

  if (ckey >= 0) return -index-1;
  return index;
}

mapkit_size_t spvector_next(spvector *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  spvector_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->key >= MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

spvector_storage *spvector_nextptr(spvector *spm, spvector_storage *pos_contents)
{
  spvector_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->key >= MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error spvector_getall(spvector *spm, spvector_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  int key;
  spvector_element *pos_array;
  spvector_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (spvector_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if (((key = pos_contents->key) >= MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error spvector_clean(spvector *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  spvector_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->key >= MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->key = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return spvector_reallocate(spm, spvector_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int spvector_compare(const void *e1, const void *e2)
{
  int key1 = ((spvector_element *)e1)->key;
  int key2 = ((spvector_element *)e2)->key;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error spvector_getall_sorted(spvector *spm, spvector_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = spvector_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), spvector_compare);

  return MAPKIT_OK;
}

mapkit_error spvector_setall(spvector *spm, spvector_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = spvector_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = spvector_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    spvector_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error spvector_removeall(spvector *spm, int *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = spvector_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  spvector_adjustcapacity(spm);

  return MAPKIT_OK;
}

void spvector_printstats(spvector *spm)
{
  fprintf(stderr, "MAPKIT: spvector statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT_spvector */



#ifdef MAPKIT__spmatrix

/*
  Implementation for _spmatrix (spmatrix_key_pair -> double)
  Default value : 0.0
  Uses a state field.
  Unassigned values always default to the default value.
*/

/* Static prototypes */

/* Return the index of key, or MAPKIT_KEYNOTFOUND */
static mapkit_size_t _spmatrix_keyindex(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash);

/* Return the index of key or -(insertion index)-1 if key not found */
static mapkit_size_t _spmatrix_insertionindex(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash);

/* Implementation */

mapkit_error _spmatrix_init(_spmatrix *spm)
{
  return _spmatrix_init_hint(spm, MAPKIT_DEFAULT_EXPECTEDUSED);
}

mapkit_error _spmatrix_init_hint(_spmatrix *spm, mapkit_size_t used)
{
#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: init\n");
#endif

  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfillfactor = 0.5;
  spm->minusedfactor = 0.2;
  spm->contents = NULL;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
  spm->defaultvalue = 0.0;
  spm->alwaysdefault = 1;

  return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, used));
}

mapkit_error _spmatrix_ensurecapacity(_spmatrix *spm, mapkit_size_t used)
{
  if (used > (spm->used + spm->maxfill - spm->fill))
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: need more capacity\n");
#endif
    return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, used));
  }
  else
    return MAPKIT_OK;
}

mapkit_error _spmatrix_adjustcapacity(_spmatrix *spm)
{
  spm->minused = (mapkit_size_t) (spm->size*spm->minusedfactor);
  spm->maxfill = (mapkit_size_t) (spm->size*spm->maxfillfactor);
  
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, spm->used));
  }
  else if (spm->fill > spm->maxfill)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
    return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, spm->used));
  }
  else
    return MAPKIT_OK;
}

void _spmatrix_free(_spmatrix *spm)
{
  free(spm->contents);
  spm->contents = NULL;
  spm->size = 0;
  spm->fill = 0;
  spm->used = 0;
  spm->maxfill = 0;
  spm->minused = 0;
#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs = spm->insertionindex_collisions = 0;
  spm->keyindexs = spm->keyindex_collisions = 0;
#endif
}

mapkit_error _spmatrix_copy(_spmatrix *to, _spmatrix *from)
{
  _spmatrix_storage *contentscopy;
  size_t size = from->size * sizeof(*from->contents);

  contentscopy = (_spmatrix_storage *)malloc(size);
  if (contentscopy == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
    
  memcpy(to, from, sizeof(_spmatrix));
  to->contents = contentscopy;
  memcpy(to->contents, from->contents, size);

  return MAPKIT_OK;
}

mapkit_size_t _spmatrix_growsize(_spmatrix *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (3*spm->minusedfactor + spm->maxfillfactor));
}

mapkit_size_t _spmatrix_shrinksize(_spmatrix *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(4.0 * used / (spm->minusedfactor + 3*spm->maxfillfactor));
}

mapkit_size_t _spmatrix_meansize(_spmatrix *spm, mapkit_size_t used)
{
  return (mapkit_size_t)(2.0 * used / (spm->minusedfactor + spm->maxfillfactor));
}

mapkit_error _spmatrix_reallocate(_spmatrix *spm, mapkit_size_t newsize)
{
  mapkit_size_t index;
  mapkit_size_t oldsize;
  _spmatrix_storage *newcontents, *oldcontents;

  /* At least one free entry */
  if (newsize <= spm->used) newsize = spm->used + 1;
  newsize = mapkit_nextprime(newsize);
  if (newsize <= spm->used)
    MAPKIT_ERROR(MAPKIT_ETOOBIG);

  newcontents = ( _spmatrix_storage *)malloc(newsize * sizeof(*spm->contents));
  if (newcontents == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);

  /* Initialize all entries to "free" */
  for (index = 0 ; index < newsize ; index ++)
    newcontents[index].state = MAPKIT_FREESLOT;

  oldcontents = spm->contents;
  oldsize = spm->size;
  spm->contents = newcontents;
  spm->size = newsize;

#ifdef MAPKIT_DEBUG
  fprintf(stderr, "MAPKIT: reallocate %ld -> %ld\n", (long)oldsize,
    (long)newsize);
#endif

  spm->maxfill = (mapkit_size_t) (newsize*spm->maxfillfactor);
  /* At least one free entry */
  if (spm->maxfill >= newsize)
    spm->maxfill = newsize - 1;
  spm->minused = (mapkit_size_t) (newsize*spm->minusedfactor);
  spm->used = 0;


  if (oldcontents != NULL)
  {
    int used = 0;
    spmatrix_key_pair key;
    double defaultvalue = spm->defaultvalue;
    int notalwaysdefault = ! spm->alwaysdefault;

    /* Copy all entries from old to new */
    for (index = 0 ; index < oldsize ; index ++)
      if (oldcontents[index].state == MAPKIT_FULLSLOT)
      {
        mapkit_size_t ins_index;
        mapkit_hash_t hash;
        _spmatrix_storage *contents;

        key = oldcontents[index].key;

        /* Fast path */
        ins_index = (hash = (mapkit_hash_t) ((key).key1 + (key).key2 * 65519)) % spm->size;
        contents = &(newcontents[ins_index]);

        if (contents->state != MAPKIT_FREESLOT)
        {
          ins_index = _spmatrix_insertionindex(spm, key, hash);
          contents = &(newcontents[ins_index]);
        }
#ifdef MAPKIT_COLLISIONS
        else
          spm->insertionindexs ++;
#endif
        if (notalwaysdefault ||
          (! ((oldcontents[index].value) == (defaultvalue))))
        {
          contents->value = oldcontents[index].value;
          contents->state = MAPKIT_FULLSLOT;
          contents->key = key;
          used ++;
        }
      }
    free(oldcontents);
    spm->used = used;
  }
  spm->fill = spm->used;

  return MAPKIT_OK;
}

double _spmatrix_value_s(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = _spmatrix_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return spm->defaultvalue;
    else
      MAPKIT_FATAL_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  return spm->contents[index].value;
}

mapkit_error _spmatrix_get_s(_spmatrix *spm, spmatrix_key_pair key,
  double *value, mapkit_hash_t hash)
{
  mapkit_size_t index;
  
  index = _spmatrix_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      *value = spm->defaultvalue;
    else
      MAPKIT_ERROR(MAPKIT_EKEYNOTFOUND);
  }

  *value = spm->contents[index].value;
  return MAPKIT_OK;
}

mapkit_error _spmatrix_set_s(_spmatrix *spm, spmatrix_key_pair key,
  double value, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = _spmatrix_insertionindex(spm, key, hash);

  if (index < 0)
    /* FULLSLOT */
    spm->contents[-index-1].value = value;
  else
  {
    _spmatrix_storage *element = &(spm->contents[index]);
    int free = element->state == MAPKIT_FREESLOT;
    element->state = MAPKIT_FULLSLOT;
    element->value = value;
    element->key = key;
    spm->used ++;
    
    if (free && ((++ spm->fill) > spm->maxfill))
    {
#ifdef MAPKIT_DEBUG
      fprintf(stderr, "MAPKIT: fill > maxfill\n");
#endif
      return _spmatrix_reallocate(spm, _spmatrix_growsize(spm, spm->used));
    }
  }
  return MAPKIT_OK;
}

double *_spmatrix_insertptr_s(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = _spmatrix_insertionindex(spm, key, hash);

  if (index < 0)
    return &(spm->contents[-index-1].value);
  else
  {
    _spmatrix_storage *element = &(spm->contents[index]);
    if (element->state == MAPKIT_FREESLOT);
    {
      /* FREESLOT */
      if (spm->fill >= spm->maxfill)
      {
        mapkit_error err;
#ifdef MAPKIT_DEBUG
        fprintf(stderr, "MAPKIT: fill => maxfill before insert\n");
#endif
        /* Must reallocate -before- inserting defaultvalue */
        err = _spmatrix_reallocate(spm, _spmatrix_growsize(spm, spm->used + 1));
        if (err)
        {
          MAPKIT_ERROR_NORET(err);
          return NULL;
        }
        
        index = _spmatrix_insertionindex(spm, key, hash);
        /* FREESLOT */
        spm->contents[index].state = MAPKIT_FULLSLOT;
        spm->contents[index].key = key;
        spm->contents[index].value = spm->defaultvalue;
        spm->used ++;
        spm->fill ++;
        return &(spm->contents[index].value);
      }
      else
        spm->fill ++;
    }

    element->state = MAPKIT_FULLSLOT;
    element->key = key;
    element->value = spm->defaultvalue;
    spm->used ++;
    return &(element->value);
  }
}

double *_spmatrix_ptr_s(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = _spmatrix_keyindex(spm, key, hash);

  if ((index >= 0) &&
    ((! spm->alwaysdefault) ||
      (! ((spm->contents[index].value) == (spm->defaultvalue)))))
    return &(spm->contents[index].value);
  else
    return NULL;
}

mapkit_error _spmatrix_remove_s(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
  mapkit_size_t index;

  index = _spmatrix_keyindex(spm, key, hash);

  if (index < 0)
  {
    if (spm->alwaysdefault)
      return MAPKIT_OK;
    else
      return MAPKIT_EKEYNOTFOUND;
  }

  spm->contents[index].state = MAPKIT_DELETEDSLOT;
  spm->used --;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return _spmatrix_reallocate(spm, _spmatrix_shrinksize(spm, spm->used));
  }

  return MAPKIT_OK;
}

mapkit_size_t _spmatrix_keyindex(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;

  signed char state;
  
  index = hash % spm->size;
  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);

#ifdef MAPKIT_COLLISIONS
  spm->keyindexs ++;
#endif
  
  while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
    && (state == MAPKIT_DELETEDSLOT ||
      (! (((spm->contents[index].key).key1 == (key).key1) && ((spm->contents[index].key).key2 == (key).key2)))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->keyindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
  }

  if (state == MAPKIT_FREESLOT) return MAPKIT_KEYNOTFOUND;
  return index;
}

mapkit_size_t _spmatrix_insertionindex(_spmatrix *spm, spmatrix_key_pair key, mapkit_hash_t hash)
{
  mapkit_size_t index, decrement;
  signed char state;

#ifdef MAPKIT_COLLISIONS
  spm->insertionindexs ++;
#endif
  
  index = hash % spm->size;

  /* Fast path (largely superfluous) */
  if ((state = spm->contents[index].state) == MAPKIT_FREESLOT) return index;
  if ((state == MAPKIT_FULLSLOT)
    && (((spm->contents[index].key).key1 == (key).key1) && ((spm->contents[index].key).key2 == (key).key2))) return -index-1;

  decrement = (hash % (spm->size-2));
  decrement += (decrement == 0);
  
  while ((state == MAPKIT_FULLSLOT)
    && (! (((spm->contents[index].key).key1 == (key).key1) && ((spm->contents[index].key).key2 == (key).key2))))
  {
#ifdef MAPKIT_COLLISIONS
    spm->insertionindex_collisions ++;
#endif
    index -= decrement;
    if (index < 0) index += spm->size;
    state = spm->contents[index].state;
  }

  if (state == MAPKIT_DELETEDSLOT)
  {
    mapkit_size_t index2 = index;
    while ((state = spm->contents[index].state) != MAPKIT_FREESLOT
      && ((state == MAPKIT_DELETEDSLOT)
        || (! (((spm->contents[index].key).key1 == (key).key1) && ((spm->contents[index].key).key2 == (key).key2)))))
    {
      index -= decrement;
      if (index < 0) index += spm->size;
      state = spm->contents[index].state;
    }
    if (state == MAPKIT_FREESLOT) return index2;
  }

  if (state == MAPKIT_FULLSLOT) return -index-1;
  return index;
}

mapkit_size_t _spmatrix_next(_spmatrix *spm, mapkit_size_t index)
{
  mapkit_size_t size = spm->size;
  _spmatrix_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;

  pos_contents = &(spm->contents[++index]);
  
  for ( ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return index;

  return -1;
}

_spmatrix_storage *_spmatrix_nextptr(_spmatrix *spm, _spmatrix_storage *pos_contents)
{
  _spmatrix_storage *end = &(spm->contents[spm->size]);
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  if (pos_contents == NULL)
    pos_contents = spm->contents;
  else
  {
    pos_contents ++;
    if (pos_contents <= spm->contents)
      return NULL;
  }

  for ( ; pos_contents < end ; pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
      return pos_contents;

  return NULL;
}

mapkit_error _spmatrix_getall(_spmatrix *spm, _spmatrix_element **array, mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  _spmatrix_element *pos_array;
  _spmatrix_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (_spmatrix_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;

  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      pos_array->key = pos_contents->key;
      pos_array->value = pos_contents->value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error _spmatrix_clean(_spmatrix *spm)
{
  mapkit_size_t index, count = 0;
  mapkit_size_t size = spm->size;
  _spmatrix_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
    
  pos_contents = spm->contents;
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && ((pos_contents->value) == (defaultvalue)))
    {
      pos_contents->state = MAPKIT_DELETEDSLOT;
      count ++;
    }

  spm->used -= count;
  if (spm->used < spm->minused)
  {
#ifdef MAPKIT_DEBUG
    fprintf(stderr, "MAPKIT: used < minused\n");
#endif
    return _spmatrix_reallocate(spm, _spmatrix_meansize(spm, spm->used));
  }

  return MAPKIT_OK;
}

int _spmatrix_compare(const void *e1, const void *e2)
{
  spmatrix_key_pair key1 = ((_spmatrix_element *)e1)->key;
  spmatrix_key_pair key2 = ((_spmatrix_element *)e2)->key;

  return (((key1).key1 < (key2).key1) ? -1 : (((key1).key1 > (key2).key1) ? 1 :
  (((key1).key2 < (key2).key2) ? -1 : (((key1).key2 == (key2).key2) ? 0 : 1))));
}

mapkit_error _spmatrix_getall_sorted(_spmatrix *spm, _spmatrix_element **array,
  mapkit_size_t *count)
{
  mapkit_error err;
  
  err = _spmatrix_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), _spmatrix_compare);

  return MAPKIT_OK;
}

mapkit_error _spmatrix_setall(_spmatrix *spm, _spmatrix_element *array, mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = _spmatrix_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = _spmatrix_set(spm, array[array_index].key, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    _spmatrix_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error _spmatrix_removeall(_spmatrix *spm, spmatrix_key_pair *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = _spmatrix_remove(spm, array[array_index]);
    if (err)
      MAPKIT_ERROR(err);
  }

  _spmatrix_adjustcapacity(spm);

  return MAPKIT_OK;
}

void _spmatrix_printstats(_spmatrix *spm)
{
  fprintf(stderr, "MAPKIT: _spmatrix statistics\n");
  fprintf(stderr, "MAPKIT: alwaysdefault = %d\n", spm->alwaysdefault);
  fprintf(stderr,
    "MAPKIT: minused = %ld, maxfill = %ld\n",
    (long)spm->minused, (long)spm->maxfill);
  fprintf(stderr,
    "MAPKIT: minusedfactor = %g, maxfillfactor = %g\n",
    spm->minusedfactor, spm->maxfillfactor);
#ifdef MAPKIT_COLLISIONS
  fprintf(stderr,
    "MAPKIT: insertionindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->insertionindexs,
    (unsigned long)spm->insertionindex_collisions);
  fprintf(stderr,
    "MAPKIT: keyindexs = %lu, collisions = %lu\n",
    (unsigned long)spm->keyindexs, (unsigned long)spm->keyindex_collisions);
#endif
}

#endif /* MAPKIT__spmatrix */



#ifdef MAPKIT_spmatrix

/* Prototypes */

/* compare the key of two spmatrix_element (for qsort) key1 then key2 */
static int spmatrix_compare1(const void *e1, const void *e2);

/* compare the key of two spmatrix_element (for qsort) key2 then key1 */
static int spmatrix_compare2(const void *e1, const void *e2);

/* Implementation */

mapkit_error spmatrix_init(spmatrix *spm)
{
  return _spmatrix_init(spm);
}

void spmatrix_free(spmatrix *spm)
{
  _spmatrix_free(spm);
}

mapkit_error spmatrix_copy(spmatrix *to, spmatrix *from)
{
  return _spmatrix_copy(to, from);
}

mapkit_error spmatrix_init_hint(spmatrix *spm,
  mapkit_size_t newsize)
{
  return _spmatrix_init_hint(spm, newsize);
}

mapkit_error spmatrix_ensurecapacity(spmatrix *spm,
  mapkit_size_t used)
{
  return _spmatrix_ensurecapacity(spm, used);
}

mapkit_error spmatrix_adjustcapacity(spmatrix *spm)
{
  return _spmatrix_adjustcapacity(spm);
}

mapkit_error spmatrix_getall(spmatrix *spm, spmatrix_element **array,
  mapkit_size_t *count)
{
  mapkit_size_t index;
  mapkit_size_t size = spm->size, vcount = 0;
  spmatrix_key_pair key;
  spmatrix_element *pos_array;
  spmatrix_storage *pos_contents;
  double defaultvalue = spm->defaultvalue;
  int notalwaysdefault = ! spm->alwaysdefault;
  
  pos_array = *array = (spmatrix_element *)malloc(sizeof(**array)*spm->used);
  if (*array == NULL)
    MAPKIT_ERROR(MAPKIT_ENOMEM);
  
  pos_contents = spm->contents;
  
  for (index = 0 ; index < size ; index ++, pos_contents ++)
    if ((pos_contents->state == MAPKIT_FULLSLOT)
      && (notalwaysdefault ||
        (! ((pos_contents->value) == (defaultvalue)))))
    {
      key = spm->contents[index].key;
      pos_array->key1 = ((key).key1);
      pos_array->key2 = ((key).key2);
      pos_array->value = spm->contents[index].value;
      pos_array ++;
      vcount ++;
    }
  *count = vcount;

  return MAPKIT_OK;
}

mapkit_error spmatrix_clean(spmatrix *spm)
{
  return _spmatrix_clean(spm);
}

int spmatrix_compare1(const void *e1, const void *e2)
{
  int key1 = ((spmatrix_element *)e1)->key1,
    key2 = ((spmatrix_element *)e2)->key1;
  int comp = ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));

  if (comp) return comp;

  key1 = ((spmatrix_element *)e1)->key2;
  key2 = ((spmatrix_element *)e2)->key2;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

int spmatrix_compare2(const void *e1, const void *e2)
{
  int key1 = ((spmatrix_element *)e1)->key2,
    key2 = ((spmatrix_element *)e2)->key2;
  int comp = ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));

  if (comp) return comp;

  key1 = ((spmatrix_element *)e1)->key1;
  key2 = ((spmatrix_element *)e2)->key1;

  return ((key1)<(key2) ? -1 : ((key1) == (key2) ? 0 : 1));
}

mapkit_error spmatrix_getall_sorted1(spmatrix *spm,
  spmatrix_element **array, mapkit_size_t *count)
{
  mapkit_error err;
  
  err = spmatrix_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), spmatrix_compare1);

  return MAPKIT_OK;
}

mapkit_error spmatrix_getall_sorted2(spmatrix *spm,
  spmatrix_element **array, mapkit_size_t *count)
{
  mapkit_error err;
  
  err = spmatrix_getall(spm, array, count);
  if (err)
    MAPKIT_ERROR(err);
  
  qsort(*array, *count, sizeof(**array), spmatrix_compare2);

  return MAPKIT_OK;
}

mapkit_error spmatrix_setall(spmatrix *spm, spmatrix_element *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;
    
  err = spmatrix_ensurecapacity(spm, spm->used + count);
  if (err)
    MAPKIT_ERROR(err);

  if (spm->alwaysdefault)
    /* Prevent shrinking */
    spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = spmatrix_set(spm, array[array_index].key1,
      array[array_index].key2, array[array_index].value);
    if (err)
      MAPKIT_ERROR(err);
  }

  if (spm->alwaysdefault)
    spmatrix_adjustcapacity(spm);

  return MAPKIT_OK;
}

mapkit_error spmatrix_removeall(spmatrix *spm, spmatrix_key *array,
  mapkit_size_t count)
{
  mapkit_size_t array_index;
  mapkit_error err;

  /* Prevent shrinking */
  spm->minused = 0;
  
  for (array_index = 0 ; array_index < count ; array_index ++)
  {
    err = spmatrix_remove(spm, array[array_index].key1, array[array_index].key2);
    if (err)
      MAPKIT_ERROR(err);
  }

  spmatrix_adjustcapacity(spm);

  return MAPKIT_OK;
}

void spmatrix_printstats(spmatrix *spm)
{
  _spmatrix_printstats(spm);
}

#endif /* MAPKIT_spmatrix */

