/* The authors of this work have released all rights to it and placed it
in the public domain under the Creative Commons CC0 1.0 waiver
(http://creativecommons.org/publicdomain/zero/1.0/).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Hash_table_(C)?oldid=19560
 */

#include"hashtbl.h"
#include<string.h>
#include<stdio.h>

static char *mystrdup(const char *s) {
    char *b;
    if (!(b = malloc(strlen(s)))) return NULL;
    strcpy(b, s);
    return b;
}

static hash_size def_hashfunc(const char *key) {
    hash_size hash = 0;

    while (*key) hash += (unsigned char) *key++;

    return hash;
}

HASHTBL *hashtbl_create(hash_size size, hash_size(*hashfunc)(const char *)) {
    HASHTBL *hashtbl;

    if (!(hashtbl = malloc(sizeof (HASHTBL)))) return NULL;

    if (!(hashtbl->nodes = calloc(size, sizeof (struct hashnode_s*)))) {
        free(hashtbl);
        return NULL;
    }

    hashtbl->size = size;

    if (hashfunc) hashtbl->hashfunc = hashfunc;
    else hashtbl->hashfunc = def_hashfunc;

    return hashtbl;
}

void hashtbl_destroy(HASHTBL *hashtbl) {
    hash_size n;
    struct hashnode_s *node, *oldnode;

    for (n = 0; n < hashtbl->size; ++n) {
        node = hashtbl->nodes[n];
        while (node) {
            free(node->key);
            oldnode = node;
            node = node->next;
            free(oldnode);
        }
    }
    free(hashtbl->nodes);
    free(hashtbl);
}

int hashtbl_insert(HASHTBL *hashtbl, const char *key, void *data) {
    struct hashnode_s *node;
    hash_size hash = hashtbl->hashfunc(key) % hashtbl->size;

    //fprintf(stderr, "hashtbl_insert() key=%s, hash=%d, data=%s\n", key, hash, (char*)data);

    node = hashtbl->nodes[hash];
    while (node) {
        if (!strcmp(node->key, key)) {
            return 0;
        }
        node = node->next;
    }

    if (!(node = malloc(sizeof (struct hashnode_s)))) return -1;
    if (!(node->key = mystrdup(key))) {
        //printf("key: %s\t data: %s\t\n", key, data);
        free(node);
        return -1;
    }
    node->data = data;
    node->next = hashtbl->nodes[hash];
    hashtbl->nodes[hash] = node;
    //printf("key: %s\t data: %s\t\n", key, data);
    return 1;
}

int hashtbl_remove(HASHTBL *hashtbl, const char *key) {
    struct hashnode_s *node, *prevnode = NULL;
    hash_size hash = hashtbl->hashfunc(key) % hashtbl->size;

    node = hashtbl->nodes[hash];
    while (node) {
        if (!strcmp(node->key, key)) {
            free(node->key);
            if (prevnode) prevnode->next = node->next;
            else hashtbl->nodes[hash] = node->next;
            free(node);
            return 0;
        }
        prevnode = node;
        node = node->next;
    }

    return -1;
}

void *hashtbl_get(HASHTBL *hashtbl, const char *key) {
    struct hashnode_s *node;
    hash_size hash = hashtbl->hashfunc(key) % hashtbl->size;

    //printf("hashtbl_get() key=%s, hash=%d\n", key, hash);

    node = hashtbl->nodes[hash];
    while (node) {
        if (!strcmp(node->key, key)) {
            //fprintf(stderr, "hashtbl_get() key=%s, data=%s\n", key, node->data);
           // printf("hashtbl_get: %s\n", node->data);
            return node->data;
        }
        node = node->next;
    }
    //printf("hashtbl_get: *NULL\n");
    return NULL;
}

int hashtbl_resize(HASHTBL *hashtbl, hash_size size) {
    HASHTBL newtbl;
    hash_size n;
    struct hashnode_s *node;

    newtbl.size = size;
    newtbl.hashfunc = hashtbl->hashfunc;

    if (!(newtbl.nodes = calloc(size, sizeof (struct hashnode_s*)))) return -1;

    for (n = 0; n < hashtbl->size; ++n) {
        for (node = hashtbl->nodes[n]; node; node = node->next) {
            hashtbl_insert(&newtbl, node->key, node->data);
            hashtbl_remove(hashtbl, node->key);

        }
    }

    free(hashtbl->nodes);
    hashtbl->size = newtbl.size;
    hashtbl->nodes = newtbl.nodes;

    return 0;
}