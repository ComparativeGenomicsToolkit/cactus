/***********************************************************************\
*
* $Source: /home/torsten/cvs/bar/bar/lists.c,v $
* $Revision: 1.2 $
* $Author: torsten $
* Contents: dynamic list functions
* Systems: all
*
\***********************************************************************/

/****************************** Includes *******************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "global.h"

#include "lists.h"

/****************** Conditional compilation switches *******************/

/***************************** Constants *******************************/

/***************************** Datatypes *******************************/

/***************************** Variables *******************************/

/****************************** Macros *********************************/

/***************************** Forwards ********************************/

/***************************** Functions *******************************/

#ifdef __cplusplus
#  extern "C" {
#endif

Node *List_newNode(ulong size)
{
  return (Node*)malloc(size);
}

Node *List_deleteNode(Node *node)
{
  Node *nextNode;

  assert(node != NULL);

  nextNode = node->next;
  free(node);

  return nextNode;
}

void List_init(void *list)
{
  assert(list != NULL);

  ((List*)list)->head  = NULL;
  ((List*)list)->tail  = NULL;
  ((List*)list)->count = 0;
}

void List_done(void                 *list,
               ListNodeFreeFunction listNodeFreeFunction,
               void                 *listNodeFreeUserData
              )
{
  assert(list != NULL);

  List_clear(list,listNodeFreeFunction,listNodeFreeUserData);
}

List *List_new(void)
{
  List *list;

  list = (List*)malloc(sizeof(List));
  if (list == NULL) return NULL;

  List_init(list);

  return list;
}

List *List_duplicate(const void           *fromList,
                     const void           *fromListFromNode,
                     const void           *fromListToNode,
                     ListNodeCopyFunction listNodeCopyFunction,
                     void                 *listNodeCopyUserData
                    )
{
  List *list;

  assert(fromList != NULL);
  assert(listNodeCopyFunction != NULL);

  list = (List*)malloc(sizeof(List));
  if (list == NULL) return NULL;

  List_init(list);
  List_copy(fromList,
            list,
            fromListFromNode,
            fromListToNode,
            NULL,
            listNodeCopyFunction,
            listNodeCopyUserData
           );

  return list;
}

void List_delete(void                 *list,
                 ListNodeFreeFunction listNodeFreeFunction,
                 void                 *listNodeFreeUserData
                )
{
  assert(list != NULL);

  List_done(list,listNodeFreeFunction,listNodeFreeUserData);\
  free(list);
}

void List_clear(void                 *list,
                ListNodeFreeFunction listNodeFreeFunction,
                void                 *listNodeFreeUserData
               )
{
  Node *node;

  assert(list != NULL);

  if (listNodeFreeFunction != NULL)
  {
    while (((List*)list)->head != NULL)
    {
      node = ((List*)list)->head;
      ((List*)list)->head = ((List*)list)->head->next;
      listNodeFreeFunction(node,listNodeFreeUserData);
    }
  }
  else
  {
    while (((List*)list)->head != NULL)
    {
      node = ((List*)list)->head;
      ((List*)list)->head = ((List*)list)->head->next;
      LIST_DELETE_NODE(node);
    }
  }
  ((List*)list)->tail  = NULL;
  ((List*)list)->count = 0;
}

void List_copy(const void           *fromList,
               void                 *toList,
               const void           *fromListFromNode,
               const void           *fromListToNode,
               void                 *toListNextNode,
               ListNodeCopyFunction listNodeCopyFunction,
               void                 *listNodeCopyUserData
              )
{
  Node *node;
  Node *newNode;
  
  assert(fromList != NULL);
  assert(toList != NULL);
  assert(listNodeCopyFunction != NULL);

  if (fromListFromNode == LIST_START) fromListFromNode = ((List*)fromList)->head;

  node = (Node*)fromListFromNode;
  while (node != fromListToNode)
  {
    newNode = listNodeCopyFunction(node,listNodeCopyUserData);
    List_insert(toList,newNode,toListNextNode);
    node = node->next;
  }
  if (node != NULL)
  {
    newNode = listNodeCopyFunction(node,listNodeCopyUserData);
    List_insert(toList,newNode,toListNextNode);
  }
}

void List_move(void *fromList,
               void *toList,
               void *fromListFromNode,
               void *fromListToNode,
               void *toListNextNode
              )
{
  Node *node;
  Node *nextNode;
  
  assert(fromList != NULL);
  assert(toList != NULL);

  if (fromListFromNode == LIST_START) fromListFromNode = ((List*)fromList)->head;

  node = (Node*)fromListFromNode;
  while (node != fromListToNode)
  {
    nextNode = node->next;
    List_remove(fromList,node);
    List_insert(toList,node,toListNextNode);
    node = nextNode;
  }
  if (node != NULL)
  {
    List_remove(fromList,node);
    List_insert(toList,node,toListNextNode);
  }
}

bool List_empty(const void *list)
{
  assert(list != NULL);
  assert(((((List*)list)->count == 0) && (((List*)list)->head == NULL) && (((List*)list)->tail == NULL)) ||
         ((((List*)list)->count > 0) && (((List*)list)->head != NULL) && (((List*)list)->tail != NULL))
        );

  return (((List*)list)->count == 0);
}

unsigned long List_count(const void *list)
{
  assert(list != NULL);
  assert(((((List*)list)->count == 0) && (((List*)list)->head == NULL) && (((List*)list)->tail == NULL)) ||
         ((((List*)list)->count > 0) && (((List*)list)->head != NULL) && (((List*)list)->tail != NULL))
        );

  return ((List*)list)->count;
}

void List_insert(void *list,
                 void *node,
                 void *nextNode
                )
{
  assert(list != NULL);

  assert(((((List*)list)->count == 0) && (((List*)list)->head == NULL) && (((List*)list)->tail == NULL)) ||
         ((((List*)list)->count > 0) && (((List*)list)->head != NULL) && (((List*)list)->tail != NULL))
        );
  assert(node != NULL);

  if      (nextNode != NULL)
  {
    ((Node*)node)->prev = ((Node*)nextNode)->prev;
    ((Node*)node)->next = ((Node*)nextNode);
    if (((Node*)nextNode)->prev != NULL) ((Node*)nextNode)->prev->next = node;
    ((Node*)nextNode)->prev = node;

    if (((List*)list)->head == nextNode) ((List*)list)->head = node;
    ((List*)list)->count++;
  }
  else if (((List*)list)->head != NULL)
  {
    ((Node*)node)->prev = ((List*)list)->tail;
    ((Node*)node)->next = NULL;

    ((List*)list)->tail->next = node;
    ((List*)list)->tail = node;
    ((List*)list)->count++;
  }
  else
  {
    ((Node*)node)->prev = NULL;
    ((Node*)node)->next = NULL;

    ((List*)list)->head  = node;
    ((List*)list)->tail  = node;
    ((List*)list)->count = 1;
  }

  assert(((((List*)list)->count == 0) && (((List*)list)->head == NULL) && (((List*)list)->tail == NULL)) ||
         ((((List*)list)->count > 0) && (((List*)list)->head != NULL) && (((List*)list)->tail != NULL))
        );
}

void List_append(void *list,
                 void *node
                )
{
  assert(list != NULL);
  assert(node != NULL);

  List_insert(list,node,NULL);
}

void *List_remove(void *list,
                  void *node
                 )
{
  void *nextNode;

  assert(list != NULL);
  assert(((List*)list)->head != NULL);
  assert(((List*)list)->tail != NULL);
  assert(((List*)list)->count > 0);
  assert((Node*)node != NULL);

  nextNode = ((Node*)node)->next;
  if (((Node*)node)->prev != NULL) ((Node*)node)->prev->next = ((Node*)node)->next;
  if (((Node*)node)->next != NULL) ((Node*)node)->next->prev = ((Node*)node)->prev;
  if ((Node*)node == ((List*)list)->head) ((List*)list)->head = ((Node*)node)->next;
  if ((Node*)node == ((List*)list)->tail) ((List*)list)->tail = ((Node*)node)->prev;
  ((List*)list)->count--;

  assert(((((List*)list)->count == 0) && (((List*)list)->head == NULL) && (((List*)list)->tail == NULL)) ||
         ((((List*)list)->count > 0) && (((List*)list)->head != NULL) && (((List*)list)->tail != NULL))
        );

  return nextNode;
}

Node *List_getFirst(void *list)
{
  Node *node;

  assert(list != NULL);

  node = ((List*)list)->head;
  if (node != NULL) List_remove(list,node);

  return node;
}

Node *List_getLast(void *list)
{
  Node *node;

  assert(list != NULL);

  node = ((List*)list)->tail;
  if (node != NULL) List_remove(list,node);

  return node;
}

const Node *List_findFirst(const void             *list,
                           ListNodeEqualsFunction listNodeEqualsFunction,
                           void                   *listNodeEqualsUserData
                          )
{
  Node *node;

  assert(list != NULL);
  assert(listNodeEqualsFunction != NULL);

  node = ((List*)list)->head;
  while ((node != NULL) && (listNodeEqualsFunction(node,listNodeEqualsUserData) != 0))
  {
    node = node->next;
  }

  return node;
}

const Node *List_findNext(const void             *list,
                          const void             *node,
                          ListNodeEqualsFunction listNodeEqualsFunction,
                          void                   *listNodeEqualsUserData
                         )
{
  assert(list != NULL);
  assert(listNodeEqualsFunction != NULL);

  UNUSED_VARIABLE(list);

  if (node != NULL)
  {
    node = (((Node*)node))->next;
    while ((node != NULL) && (listNodeEqualsFunction(node,listNodeEqualsUserData) != 0))
    {
      node = (((Node*)node))->next;
    }
  }

  return node;
}

#if 0
void pp(void *list)
{
  void *node;

printf("---\n");
  node = ((List*)list)->head;
  while (node != NULL)
  {
printf("%p\n",node);
node = ((Node*)node)->next;
  }
}
#endif /* 0 */

void List_sort(void                    *list,
               ListNodeCompareFunction listNodeCompareFunction,
               void                    *listNodeCompareUserData
              )
{
  List  sortedList;
  void  *node1,*node2;
  ulong n;
  bool  mergedFlag;
  ulong i;
  ulong n1,n2;
  void  *node;

  assert(list != NULL);
  assert(listNodeCompareFunction != NULL);

//pp(list);

  /* sort list with merge-sort */
  n = 1;
  do
  {
    sortedList.head = NULL;
    sortedList.tail = NULL;

    mergedFlag = FALSE;
    node1 = ((List*)list)->head;
    while (node1 != NULL)
    {
      /* find start of sub-list 2 */
      node2 = node1;
      for (i = 0; (i < n) && (node2 != NULL); i++)
      {
        node2 = ((Node*)node2)->next;
      }

      /* merge */
      n1 = n;
      n2 = n;
      while (((n1 > 0) && (node1 != NULL)) || ((n2 > 0) && (node2 != NULL)))
      {
        /* select next node to add to sorted list */
        if      ((n1 == 0) || (node1 == NULL))
        {
          /* sub-list 1 is empty -> select node from sub-list 2 */
          node = node2; node2 = ((Node*)node2)->next; n2--;
        }
        else if ((n2 == 0) || (node2 == NULL))
        {
          /* sub-list 2 is empty -> select node from sub-list 1 */
          node = node1; node1 = ((Node*)node1)->next; n1--;
        }
        else
        {
          /* compare nodess from sub-list 1, 2 */
          if (listNodeCompareFunction(node1,node2,listNodeCompareUserData) < 0)
          {
            /* node1 < node2 -> select node1 */
            node = node1; node1 = ((Node*)node1)->next; n1--;
          }
          else
          {
            /* node1 >= node2 -> select node2 */
            node = node2; node2 = ((Node*)node2)->next; n2--;
          }
          mergedFlag = TRUE;    
        }

        /* add to list */
        ((Node*)node)->prev = ((List*)list)->tail;
        ((Node*)node)->next = NULL;
        if (sortedList.head != NULL)
        {
          sortedList.tail->next = node;
          sortedList.tail = node;
        }
        else
        {
          sortedList.head = node;
          sortedList.tail = node;
        }
      }
//pp(&sortedList);

      /* next sub-lists */
      node1 = node2;
    }

    /* next sub-list size */
    ((List*)list)->head = sortedList.head;
    ((List*)list)->tail = sortedList.tail;
    n *= 2;
  }
  while (mergedFlag);
}

#ifdef __cplusplus
  }
#endif

/* end of file */
