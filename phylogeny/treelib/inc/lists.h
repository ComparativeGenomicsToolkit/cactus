/***********************************************************************\
*
* $Source: /home/torsten/cvs/bar/bar/lists.h,v $
* $Revision: 1.2 $
* $Author: torsten $
* Contents: dynamic list functions
* Systems: all
*
\***********************************************************************/

#ifndef __LISTS__
#define __LISTS__

/****************************** Includes *******************************/
#include <stdlib.h>

#include "global.h"

/****************** Conditional compilation switches *******************/

/***************************** Constants *******************************/

#define LIST_START NULL
#define LIST_END   NULL

/***************************** Datatypes *******************************/

#define LIST_NODE_HEADER(type) \
  type *prev; \
  type *next

#define LIST_HEADER(type) \
  type          *head; \
  type          *tail; \
  unsigned long count

typedef struct Node
{
  LIST_NODE_HEADER(struct Node);
} Node;

typedef struct
{
  LIST_HEADER(Node);
} List;

/* delete list node function */
typedef void(*ListNodeFreeFunction)(void *node, void *userData);

/* copy list node function */
typedef void*(*ListNodeCopyFunction)(const void *node, void *userData);

/* list node equals function */
typedef int(*ListNodeEqualsFunction)(const void *node, void *userData);

/* compare list nodes function */
typedef int(*ListNodeCompareFunction)(const void *node1, const void *node2, void *userData);

/***************************** Variables *******************************/

/****************************** Macros *********************************/

#define LIST_STATIC_INIT {NULL,NULL}

#define LIST_NEW_NODE(Type) (Type*)List_newNode(sizeof(Type))
#define LIST_DELETE_NODE(node) List_deleteNode((Node*)node)

#define LIST_DEFINE(type,define) \
  typedef struct { define; } type; \
  typedef struct type ## Node\
  { \
    LIST_NODE_HEADER(struct type ## Node); \
    define; \
  } type ## Node; \
  typedef struct \
  { \
    LIST_HEADER(type ## Node); \
  } type ## List

/***************************** Forwards ********************************/

/***************************** Functions *******************************/

#ifdef __cplusplus
  extern "C" {
#endif

/***********************************************************************\
* Name   : List_newNode
* Purpose: allocate new list node
* Input  : size - size of node
* Output : -
* Return : node or NULL if insufficient memory
* Notes  : -
\***********************************************************************/

Node *List_newNode(ulong size);

/***********************************************************************\
* Name   : List_deleteNode
* Purpose: delete list node
* Input  : node - list node
* Output : -
* Return : next node in list or NULL
* Notes  : -
\***********************************************************************/

Node *List_deleteNode(Node *node);

/***********************************************************************\
* Name   : List_init
* Purpose: initialise list
* Input  : list - list to initialize
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_init(void *list);

/***********************************************************************\
* Name   : List_done
* Purpose: free all nodes
* Input  : list                 - list to free
*          listNodeFreeFunction - free function for single node or NULL
*          listNodeFreeUserData - user data for free function
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_done(void                 *list,
               ListNodeFreeFunction listNodeFreeFunction,
               void                 *listNodeFreeUserData
              );

/***********************************************************************\
* Name   : List_new
* Purpose: allocate new list
* Input  : -
* Output : -
* Return : list or NULL on insufficient memory
* Notes  : -
\***********************************************************************/

List *List_new(void);

/***********************************************************************\
* Name   : List_duplicate
* Purpose: duplicate list
* Input  : fromList                        - from list
*          fromListFromNode,fromListToNode - from/to node (could be
*                                            NULL)
*          listNodeCopyFunction            - node copy function
*          listNodeCopyUserData            - node copy user data
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

List *List_duplicate(const void           *fromList,
                     const void           *fromListFromNode,
                     const void           *fromListToNode,
                     ListNodeCopyFunction listNodeCopyFunction,
                     void                 *listNodeCopyUserData
                    );

/***********************************************************************\
* Name   : List_delete
* Purpose: free all nodes and delete list
* Input  : list                 - list to free
*          listNodeFreeFunction - free function for single node or NULL
*          listNodeFreeUserData - user data for free function
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_delete(void                 *list,
                 ListNodeFreeFunction listNodeFreeFunction,
                 void                 *listNodeFreeUserData
                );

/***********************************************************************\
* Name   : List_clear
* Purpose: free all nodes in list
* Input  : list                 - list
*          listNodeFreeFunction - free function for single node or NULL
*          listNodeFreeUserData - user data for free function
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_clear(void                 *list,
                ListNodeFreeFunction listNodeFreeFunction,
                void                 *listNodeFreeUserData
               );

/***********************************************************************\
* Name   : List_copy
* Purpose: copy contents of list
* Input  : fromList                        - from list
*          toList                          - to list
*          fromListFromNode,fromListToNode - from/to node (could be
*                                            NULL)
*          toListNextNode                  - insert node before nextNode
*                                            (could be NULL)
*          listNodeCopyFunction            - node copy function
*          listNodeCopyUserData            - node copy user data
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_copy(const void           *fromList,
               void                 *toList,
               const void           *fromListFromNode,
               const void           *fromListToNode,
               void                 *toListNextNode,
               ListNodeCopyFunction listNodeCopyFunction,
               void                 *listNodeCopyUserData
              );

/***********************************************************************\
* Name   : List_move
* Purpose: move contents of list
* Input  : fromList                        - from list
*          toList                          - to list
*          fromListFromNode,fromListToNode - from/to node (could be
*                                            NULL)
*          toListNextNode                  - insert node before nextNode
*                                            (could be NULL)
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_move(void *fromList,
               void *toList,
               void *fromListFromNode,
               void *fromListToNode,
               void *toListNextNode
              );

/***********************************************************************\
* Name   : List_empty
* Purpose: check if list is empty
* Input  : list - list
* Output : -
* Return : TRUE if list is empty, FALSE otherwise
* Notes  : -
\***********************************************************************/

bool List_empty(const void *list);

/***********************************************************************\
* Name   : List_count
* Purpose: get number of elements in list
* Input  : list - list
* Output : -
* Return : number of elements
* Notes  : -
\***********************************************************************/

unsigned long List_count(const void *list);

/***********************************************************************\
* Name   : List_ins
* Purpose: insert node into list
* Input  : list     - list
*          node     - node to insert
*          nextNode - insert node before nextNode (could be NULL)
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_insert(void *list,
                 void *node,
                 void *nextNode
                );

/***********************************************************************\
* Name   : List_append
* Purpose: append node to end of list
* Input  : list - list
*          node - node to add
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void List_append(void *list,
                 void *node
                );

/***********************************************************************\
* Name   : List_remove
* Purpose: remove node from list
* Input  : list - list
*          node - node to remove
* Output : -
* Return : next node in list or NULL
* Notes  : -
\***********************************************************************/

void *List_remove(void *list,
                  void *node
                 );

/***********************************************************************\
* Name   : List_getFirst
* Purpose: remove first node from list
* Input  : list - list
* Output : -
* Return : removed node or NULL if list is empty
* Notes  : -
\***********************************************************************/

Node *List_getFirst(void *list);

/***********************************************************************\
* Name   : List_getLast
* Purpose: remove last node from list
* Input  : list - list
* Output : -
* Return : removed node or NULL if list is empty
* Notes  : -
\***********************************************************************/

Node *List_getLast(void *list);

/***********************************************************************\
* Name   : List_findFirst
* Purpose: find node in list
* Input  : list                   - list
*          listNodeEqualsFunction - equals function
*          listNodeEqualsUserData - user data for equals function
* Output : -
* Return : node or NULL if not found
* Notes  : -
\***********************************************************************/

const Node *List_findFirst(const void             *list,
                           ListNodeEqualsFunction listNodeEqualsFunction,
                           void                   *listNodeEqualsUserData
                          );

/***********************************************************************\
* Name   : List_findNext
* Purpose: find next node in list
* Input  : list                    - list
*          node                    - previous found node
*          listNodeEqualsFunction - equals function
*          listNodeEqualsUserData - user data for equals function
* Output : -
* Return : next node or NULL if no next node found
* Notes  : -
\***********************************************************************/

const Node *List_findNext(const void             *list,
                          const void             *node,
                          ListNodeEqualsFunction listNodeEqualsFunction,
                          void                   *listNodeEqualsUserData
                         );

/***********************************************************************\
* Name   : List_sort
* Purpose: sort list
* Input  : list                - list
*          listNodeCompareFunction - compare function
*          listNodeCmpUserData - user data for compare function
* Output : -
* Return : -
* Notes  : use temporary O(n) memory
\***********************************************************************/

void List_sort(void                    *list,
               ListNodeCompareFunction listNodeCompareFunction,
               void                    *listNodeCompareUserData
              );

#ifdef __cplusplus
  }
#endif

#endif /* __LISTS__ */

/* end of file */
