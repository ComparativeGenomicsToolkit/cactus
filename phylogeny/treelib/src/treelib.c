#include "treelib.h"

float
jcdist(char *seqA, char *seqB)
{
  float p = 0.0f;

  unsigned int seqA_len = strlen(seqA);
  unsigned int seqB_len = strlen(seqB);

  assert(seqA_len == seqB_len);

  unsigned int i = 0;
  char seqA_base;
  char seqB_base;

  for (i=0; i<seqA_len; i++) {
    seqA_base = toupper(seqA[i]);
    seqB_base = toupper(seqB[i]);

    if (seqA_base == 'N' || seqB_base == 'N') {
      p += 0.75f;
    } else if (seqA_base != seqB_base) {
      p += 1.0f;
    }
  }

  p /= (float) seqA_len;

  float dist = 10.0f;
  if (p < 0.75f) {
    dist = -0.75f * log(1.0f - (4.0f/3.0f) * p);
  }

  if (dist == -0.0f) {
    dist = 0.0f;
  }

  return dist;
}

struct DistanceMatrix *
mfa2dist (char **aln, unsigned int num)
{
  struct DistanceMatrix *mat;
  mat = empty_DistanceMatrix(num);

  unsigned int i, j;

  for (i=0; i<num; i++) {
    for (j=i+1; j<num; j++) {
      mat->data[i][j] = jcdist(aln[i], aln[j]);
      mat->data[j][i] = jcdist(aln[i], aln[j]);
    }
  }

  return mat;
}

void
aln_init (struct Alignment **aln_loc, unsigned int size)
{
  *aln_loc = (struct Alignment *) malloc_util( sizeof(struct Alignment ));

  (*aln_loc)->numseqs = size;
  (*aln_loc)->seqs = (struct Sequence **) malloc_util( size * sizeof( struct Sequence *));
  (*aln_loc)->length = 0;

  unsigned int i = 0;
  char *buffer = NULL;
  for (i=0; i<size; i++) {
    (*aln_loc)->seqs[i] = empty_Sequence();
    
    asprintf(&buffer, "%d", i);
    (*aln_loc)->seqs[i]->name = buffer;
  }

  return;
}

void
distance2text (String s, char *species, double dist) {
  char distbuffer[32];

  String_appendCString(s, species);
  String_appendChar(s, ':');

  sprintf(distbuffer, "%g", dist);
  String_appendCString(s, distbuffer);
}

void
tree2string_node(struct Tnode *node, String s)
{

  if (node != NULL) {
    if (node->left == NULL && node->right == NULL) {
      if (node->clust->clustersize == 0 || node->clust->members == NULL) {
      } else if (node->clust->clustersize == 1) {
        distance2text(s, node->clust->members[0]->name, node->distance);
      } else {
        unsigned int i;
        for (i=0; i<node->clust->clustersize-1; i++) {
          String_appendChar(s, '(');
          distance2text(s, node->clust->members[i]->name, 0.0);
        }
        distance2text(s, node->clust->members[i]->name, 0.0);
        String_appendChar(s, ')');
        for (i=0; i<node->clust->clustersize-2; i++) {
          String_appendChar(s, ':');
          distance2text(s, "", 0.0);
          String_appendChar(s, ')');
        }

        String_appendChar(s, ':');
        distance2text(s, "", node->distance);
      }
    } else if (node->left != NULL && node->right != NULL) {
      String_appendChar(s, '(');
      tree2string_node(node->left, s);
      String_appendChar(s, ',');
      tree2string_node(node->right, s);
      String_appendChar(s, ')');
      distance2text(s, "", node->distance);
    }
  }
}

char *
tree2string (struct Tree *tree) {
  String s;
  s = String_new();

  if (tree != NULL) {
    if (tree->child[0] != NULL) {
       if (tree->child[1] == NULL) {
         if (tree->child[0]->left != NULL && tree->child[0]->right != NULL) {
           String_appendChar(s, '(');
           tree2string_node(tree->child[0]->left, s);
           String_appendChar(s, ',');
           tree2string_node(tree->child[0]->right, s);
           String_appendCString(s, ");");           

         } else {
           if (tree->child[0]->clust->clustersize == 1) {
             String_appendCString(s, "0;");
           } else {
             unsigned  int i;
             for (i=0; i<tree->child[0]->clust->clustersize-1; i++) {
               String_appendChar(s, '(');
               distance2text(s, tree->child[0]->clust->members[i]->name, 0.0);
             }

             distance2text(s, tree->child[0]->clust->members[i]->name, 0.0);

             for (i=0; i<tree->child[0]->clust->clustersize-2; i++) {
               String_appendCString(s, "):0.0)");
             }

             String_appendCString(s, ");");
           }
         }
       }
    }
  }

  unsigned int strlen = String_length(s);
  char *buffer = NULL;
  buffer = malloc(sizeof(char) * (strlen+1));
  strncpy(buffer, String_cString(s), strlen);
  String_delete(s);

  return buffer;
}

char *
msa2tree (char **mfa, unsigned int num) {
  struct DistanceMatrix *mat = NULL;
  struct Alignment *aln = NULL;
  struct ClusterGroup *group = NULL;
  struct Tree *njTree = NULL;

  char *treestring = NULL;
  assert(num >= 0);
  if (num == 0) {
    return NULL;
  } else if (num == 1) {
    asprintf(&treestring, "%d;", 0);
    return treestring;
  } else if (num == 2) {
    float dist = jcdist(mfa[0], mfa[1]);
    dist *= 0.5;
    asprintf(&treestring, "(%d:%g, %d:%g)", 0, dist, 1, dist);
    return treestring;
  }

  mat = mfa2dist(mfa, num);
  aln_init(&aln, num);
  group = alignment_to_ClusterGroup(aln, 0);
  group->matrix = mat;

  njTree = neighbour_joining_buildtree(group, 0);

  struct Tnode *tmpNode = NULL;
  tmpNode = njTree->child[2];
  njTree->child[2] = NULL;

  struct Tnode *leftNode = NULL;
  leftNode = (struct Tnode *) malloc_util(sizeof(struct Tnode));
  leftNode->left = njTree->child[0];
  leftNode->right = njTree->child[1];
  leftNode->distance = 0.0;

  struct Tnode *rootNode = NULL;
  rootNode = (struct Tnode *) malloc_util(sizeof(struct Tnode));
  rootNode->left = leftNode;
  rootNode->right = tmpNode;

  njTree->child[0] = rootNode;
  njTree->child[1] = NULL;

  treestring = tree2string(njTree);

  aln = free_Alignment(aln);
  group = free_ClusterGroup(group);
  njTree = free_Tree(njTree);

  return treestring;
}
