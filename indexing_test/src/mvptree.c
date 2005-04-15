#include <time.h>
#include "mvptree.h"
#include "smatrix.h"
#include "ufragdb.h"

/****************************************************************************
The random function randomly selects an integer from the interval [low,high].
****************************************************************************/
static 
int rand_num(int low, int high)
{
	int result=lrand48(), interval=high-low+1;

	result=(result%interval)+low;
	return result;
}


/***********************************************************************************
The adjust function is used to heapify the arrays distlist and oidlist to be sorted.
***********************************************************************************/
static
void adjust(int *oidlist, DIST_TYPE *distlist, int first, int i, int n)
{
	int j, int_temp;
	DIST_TYPE dist, float_temp;

	first--;
	int_temp=oidlist[i+first]; 	
	float_temp=distlist[i+first];
	dist=distlist[i+first];
	j=2*i;
	while(j<=n) {
		if((j<n)&&((distlist[j+first]<distlist[j+first+1]))) j++;
		if(dist>=distlist[j+first]) break;
		distlist[(j/2)+first]=distlist[j+first];
		oidlist[(j/2)+first]=oidlist[j+first];
		j*=2;
	}
	distlist[(j/2)+first]=float_temp;
	oidlist[(j/2)+first]=int_temp;
}

/*******************************************************************
The sort function sorts  the objects with respect to their distances
to a vantage point given in distlist.  Heapsort algorithm is used.
In particular, the points referred in oidlist[first] to
oidlist[last] are sorted with respect to the distances given in
distlist[first} to distlist[last]. 
The indices first and last are used inclusively.
*******************************************************************/
static
void sort(int *oidlist, DIST_TYPE *distlist, int first, int last)
{
	int size, i, int_temp;
	DIST_TYPE float_temp;

	size=last-first+1;
	for(i=size/2; i>0; i--) adjust(oidlist, distlist, first, i, size);
	for(i=size-1; i>0; i--) {
		int_temp=oidlist[first+i];
		float_temp=distlist[first+i];
		oidlist[first+i]=oidlist[first];
		distlist[first+i]=distlist[first];
		oidlist[first]=int_temp;
		distlist[first]=float_temp;
		adjust(oidlist, distlist, first, 1, i);
 	}
}

static
void insert_hit(MVP_TREE *MVPT, int obj, DIST_TYPE d)
{
  if (MVPT->hits_used >= MVPT->hits_avail) {
    MVPT->hits_avail += MAX_HITS;
    MVPT->hits = reallocec(MVPT->hits, MVPT->hits_avail * sizeof(int)); 
    MVPT->hit_dists = reallocec(MVPT->hit_dists, 
				MVPT->hits_avail * sizeof(DIST_TYPE)); 
  }
  MVPT->hits[MVPT->hits_used] = obj;
  MVPT->hit_dists[MVPT->hits_used] = d;
  MVPT->hits_used++;
  MVPT->seqs_hit++;
}


/*********************************************************/
/* The function below computes the distance between two data 
   points (actually a data point and a vantage point).
   It is used in construction of the mvp-tree. */
static
DIST_TYPE distance(MVP_TREE *MVPT, int obj_1, int obj_2)
{
  BIOSEQ s1; 
  BIOSEQ s2; 
  MVPT->data_pt(MVPT->db, obj_2, &s1);
  MVPT->data_pt(MVPT->db, obj_1, &s2);
  MVPT->num_calc_con++;
  return (DIST_TYPE) MVPT->metric(MVPT->matrix, &s1, &s2); 
}

/*********************************************************/
/* The function below computes the distance between a data
   point and the query point. It is used in the search 
   function of the mvp-tree. */

static
DIST_TYPE dist(MVP_TREE *MVPT, BIOSEQ *q, int obj)
{
  BIOSEQ s;
  MVPT->data_pt(MVPT->db, obj, &s);
  return (DIST_TYPE) MVPT->metric(MVPT->matrix, q, &s);
}


/*********************************************************/
/* The function below returns 1 if the interval [x1,y1] 
   intersects with the interval [x2,y2]. */
static 
int intersect(DIST_TYPE x1, DIST_TYPE x2, DIST_TYPE yy1, DIST_TYPE yy2)
{
  if((x1<=yy2)&&(yy1<=x2)) return 1;
  else return 0;
}


static
int alloc_node(MVP_TREE *MVPT)
{
  int new_node;
  if (MVPT->nodes_used >= MVPT->nodes_avail) {
    MVPT->nodes_avail += NODES_INCR;
    MVPT->nodes = reallocec(MVPT->nodes, MVPT->nodes_avail * sizeof(NODE)); 
  }
  new_node = MVPT->nodes_used++;
  return new_node;
}




/**********************************************************************
The function create_vp-tree creates a vp-tree for the objects whose ids 
are listed in oidlist[] between the FIRST and LAST (inclusive) entries. 
It returns the pointer to the root of the tree. 

oidlist keeps the array indexes of data points.
Sorting with respect to the distance from the vantage points
are done on this array. The array distlist keeps the distances
from the vantage points and is only used at the time of creation. 
**********************************************************************/


static
int create_mvp_tree(MVP_TREE *MVPT, int *oidlist, DIST_TYPE *distlist, 
		    int level, int FIRST, int LAST)
{
	NODE *newnode; 
	int newnode_id;
	int vp1, vp2, temp, i, jj1, j2, k, size, first, last; 

	/* create a new node */ 
	if (FIRST == 62624)
	  printf("create: level=%d, FIRST=%d, LAST=%d\n", level, FIRST, LAST) ;

	newnode_id = alloc_node(MVPT);
	newnode = MVPT->nodes + newnode_id;
	
	/* if there are at most K+2 (K data elements and 2 vantage
	points) elements, create a leaf node. */ 
	if(LAST-FIRST <= MVPT_K+1) {

	  if (FIRST == 62624)
	    printf("I'm here 01\n");

	  /* LAST and FIRST are inclusive array bounds, that is why we
	     compare their difference to K+1 */ 
	  newnode->node_type=LEAF;
	  /* The size includes the second vantage point. In
	     other words it can at most be K+1 */
	  newnode->size=LAST-FIRST; 
	  if(LAST==FIRST) {
	    newnode->vp1=oidlist[FIRST]; 
	    return(newnode_id); 
	  } 
	  /* pick the first vantage point randomly */
	  temp=rand_num(FIRST, LAST); 
	  vp1=oidlist[temp]; 
	  /* swap the vantage point with the first element */
	  if (FIRST == 62624) {
	    printf("I'm here 011\n");
	    printf("temp=%d, FIRST = %d\n", temp, FIRST);
	  }
	  oidlist[temp]=oidlist[FIRST]; 
	  oidlist[FIRST]=vp1;
	  FIRST++; 

	  if (FIRST == 62624)
	    printf("I'm here 0111\n");

	  /* compute distances from the first vantage point */ 
	  for(i=FIRST; i<=LAST; i++) 
	    distlist[i] = distance(MVPT, vp1, oidlist[i]); 
	  
	  if (FIRST == 62624)
	    printf("I'm here 012\n");

	  sort(oidlist, distlist, FIRST,LAST); 


	  if (FIRST == 62624)
	    printf("I'm here 02\n");

	  /* pick the second vantage point to be the farthest point from the first one */ 
	  vp2=oidlist[LAST];
	  oidlist[LAST]=oidlist[FIRST]; 
	  oidlist[FIRST]=vp2;
	  newnode->vp1_to_vp2=distlist[LAST];
	  distlist[LAST]=distlist[FIRST]; 
	  FIRST++;
	  newnode->vp1=vp1; 
	  newnode->vp2=vp2; 

	  for(i=FIRST; i<=LAST; i++) {
	    newnode->CD.leaf.datapoint[i-FIRST]=oidlist[i];
	    newnode->CD.leaf.dist1[i-FIRST]=distlist[i]; 
	  }

	  if (FIRST == 62624)
	    printf("I'm here 04\n");

	  for(i=FIRST; i<=LAST; i++) { 
	    distlist[i] = distance(MVPT, vp2, oidlist[i]);
	    newnode->CD.leaf.dist2[i-FIRST]=distlist[i]; 
	  }
	  return(newnode_id); 
	} 
	/* If there are more than K+2 elements, we should create an internal node. */ 
	/* pick the first vantage point randomly */
	temp=rand_num(FIRST,LAST); 
	vp1=oidlist[temp];
	oidlist[temp]=oidlist[FIRST]; 
	oidlist[FIRST]=vp1;
	FIRST++; 
	/* the parameter level keeps track of the number of ancestor
	   vantage point above this node (the one that is being
	   created).  If it is less than MVPT->num_ref (parameter p in the
	   paper), we save the  distance from that vantage point */ 
	if(level > MVPT->num_ref)
	  for(i=FIRST; i<=LAST; i++) 
	    distlist[i] = distance(MVPT, vp1, oidlist[i]); 
	else 
	  for(i=FIRST; i<=LAST; i++) {
	    distlist[i] = distance(MVPT, vp1, oidlist[i]); 
	    MVPT->refdist[oidlist[i]][level-1] = distlist[i]; 
	  } 
	level++;
	sort(oidlist, distlist, FIRST, LAST); 
	/* Now the distances from vp1 are sorted and kept in distlist[].
	   We have to find the cutoff points next. */
	newnode->CD.internal.cutoff1[MVPT_M] = distlist[LAST];
	newnode->node_type=INTERNAL; newnode->vp1=vp1;
	newnode->size = MVPT_M * MVPT_N; 
	/* select the second vantage point from one of the farthest
	   points to the first one. */ 
	temp=rand_num(LAST-MVPT_N+1, LAST); 
	vp2=oidlist[temp];
	newnode->vp2=vp2; 
	newnode->vp1_to_vp2=distlist[temp];
	for(i=temp; i<LAST; i++) { 
	  oidlist[i]=oidlist[i+1];
	  distlist[i]=distlist[i+1]; 
	} 
	oidlist[LAST]=vp2;
	LAST--; 
	size=LAST-FIRST+1; 
	jj1=size/MVPT_M; 
	/* Create M partitions with respect to vp1 */ 
	for(i=0; i < MVPT_M; i++) {
	  first=FIRST+(i*jj1);
	  newnode->CD.internal.cutoff1[i]=distlist[first];
	  if(i == MVPT_M-1) 
	    size=LAST-FIRST-(i*jj1)+1; 
	  else size=jj1;
	  last=first+size-1; 
	  if(level > MVPT->num_ref) 
	    for(k=first; k<=last; k++) 
	      distlist[k] = distance(MVPT, vp2, oidlist[k]);
	  else for(k=first; k<=last; k++) { 
	    distlist[k] = distance(MVPT, vp2, oidlist[k]);
	    MVPT->refdist[oidlist[k]][level-1] = distlist[k]; 
	  }
	  sort(oidlist, distlist, first, last);
	  newnode->CD.internal.cutoff2[i][MVPT_N]=distlist[last];
	  j2=size/MVPT_N; 
	  /* For each of the M partitions, create N partitions with
	     respect to vp2 */ 
	  for(k=0; k < MVPT_N; k++) {
	    newnode->CD.internal.cutoff2[i][k]=distlist[first+j2*k];
	    if (k== MVPT_N-1)
	      newnode->CD.internal.child[i][k] = 
		create_mvp_tree(MVPT, oidlist, distlist, level+1, 
				first+j2*k, last); 
	    else
	      newnode->CD.internal.child[i][k] = 
		create_mvp_tree(MVPT, oidlist, distlist, level+1, 
				first+j2*k, first+j2*k-1+j2); 
	  } 
	}
	return(newnode_id); 
}




/*********************************************************/
/* The search function looks for data points that are within
   r of the query object. */
static 
void search(MVP_TREE *MVPT, int level, int node_id, BIOSEQ *queryobj, 
	    DIST_TYPE r)
{
  int i, j, k;
  DIST_TYPE dist_from_vp1, dist_from_vp2, d;
  char calc_dist;
  NODE *node = MVPT->nodes + node_id;

  /* compute the distance between the query point and the first
     vantage point of the node */ 
  dist_from_vp1 = dist(MVPT, queryobj, node->vp1);
  MVPT->seqs_visited++;

  if(level <= MVPT->num_ref) 
    MVPT->SPdist[level-1]=dist_from_vp1;
  level++;

  if(dist_from_vp1 <= r) 
    insert_hit(MVPT, node->vp1, dist_from_vp1); // vp1 is in the answer set

  if(node->node_type==LEAF) {
    if(node->size==0) return;
    if(node->size==1) { 
      if((node->vp1_to_vp2 > dist_from_vp1 + r) || 
	 (node->vp1_to_vp2 < dist_from_vp1 - r))   
	return;
      /* the loop below makes use of the stored distance computations
	 from a data point to the first MVPT->num_ref top level ancestor
	 vantage points. */ 
      for (k=0; k < MVPT->num_ref; k++)
	if((MVPT->refdist[node->vp2][k] > MVPT->SPdist[k] + r) ||
	   (MVPT->refdist[node->vp2][k] < MVPT->SPdist[k] - r)) 
	  return;
    }

    /* compute the distance from the query object to the second vantage point */
    dist_from_vp2 = dist(MVPT, queryobj, node->vp2);
    MVPT->seqs_visited++;

    if(dist_from_vp2 <= r) 
      insert_hit(MVPT, node->vp2, dist_from_vp2); // vp2 is in the answer set

    /* check for each data point in the leaf */
    for (i=0; i < node->size-1; i++) 
      if((node->CD.leaf.dist1[i] <= dist_from_vp1 + r) && 
	 (node->CD.leaf.dist1[i] >= dist_from_vp1 - r) &&
	 (node->CD.leaf.dist2[i] <= dist_from_vp2 + r) &&
	 (node->CD.leaf.dist2[i] >= dist_from_vp2 - r)) {
	calc_dist=0;
	for (k=0; k < MVPT->num_ref; k++)
	  if((MVPT->refdist[node->CD.leaf.datapoint[i]][k] > MVPT->SPdist[k] + r) || 
	     (MVPT->refdist[node->CD.leaf.datapoint[i]][k] < MVPT->SPdist[k] - r)) {
	    calc_dist=1;
	    break;
	  }
	if(calc_dist) continue;
	d = dist(MVPT, queryobj, node->CD.leaf.datapoint[i]);
	MVPT->seqs_visited++;
	if(d <= r) 
	  insert_hit(MVPT, node->CD.leaf.datapoint[i], d);
      }
    return;
  }

  /* else if node->node_type==INTERNAL */
  dist_from_vp2 = dist(MVPT, queryobj, node->vp2);
  MVPT->seqs_visited++;

  if(level <= MVPT->num_ref) 
    MVPT->SPdist[level-1] = dist_from_vp2;
  level++;

  if(dist_from_vp2 <= r) 
    insert_hit(MVPT, node->vp2, dist_from_vp2); // vp2 is in the answer set

  /* Continue to search following the children of the current node */
  for(i=0; i < MVPT_M; i++) {
    if(intersect(node->CD.internal.cutoff1[i], 
		 node->CD.internal.cutoff1[i+1], 
		 dist_from_vp1 - r, 
		 dist_from_vp1 + r)) {  
      for(j=0; j < MVPT_N; j++) {
	MVPT->nodes_visited++;
	if (intersect(node->CD.internal.cutoff2[i][j], 
		      node->CD.internal.cutoff2[i][j+1], 
		      dist_from_vp2 - r,
		      dist_from_vp2 + r)) {
	  MVPT->nodes_hit++;
	  search(MVPT, level, node->CD.internal.child[i][j], queryobj, r);
	}
      }
    }
  }
}



MVP_TREE *MVP_TREE_init(void *db, void *matrix,  dpt_func *data_pt, 
			mtr_func *metric, int num_ref, 
			DIST_TYPE **refdist)
{
  MVP_TREE *MVPT = callocec(1, sizeof(MVP_TREE));
  
  MVPT->db = db;
  MVPT->matrix = matrix;
  MVPT->data_pt = data_pt;
  MVPT->metric = metric;
  MVPT->num_ref = num_ref;
  MVPT->refdist = refdist;

  return MVPT;
}


void MVP_TREE_create(MVP_TREE *MVPT, int *oid, int oid_size)
{
  DIST_TYPE *distlist;

  if (MVPT == NULL) return;

  MVPT->no_pts = oid_size;

  if (MVPT->num_ref > 0) {
    MVPT->SPdist = mallocec(MVPT->num_ref * sizeof(DIST_TYPE));
  }

  srand48(time(NULL));
  distlist = mallocec(oid_size * sizeof(DIST_TYPE));

  MVPT->nodes = mallocec(NODES_INCR * sizeof(NODE));
  /* The root is always 0 */
  printf("Before create_mvp_tree\n");
  create_mvp_tree(MVPT, oid, distlist, 1, 0, oid_size-1);
  printf("After create_mvp_tree\n");
  free(distlist);
}

void MVP_TREE_destroy(MVP_TREE *MVPT)
{
  if (MVPT->num_ref > 0) {
    free(MVPT->SPdist);
  }

  free(MVPT->hits);
  free(MVPT->hit_dists);
  free(MVPT->nodes);
  free(MVPT);
}

int MVP_TREE_write(MVP_TREE *MVPT, FILE *fp)
{
  int size = 0;
  size += sizeof(int) * fwrite(&MVPT->no_pts, sizeof(int), 1, fp);
  size += sizeof(ULINT) * fwrite(&MVPT->num_calc_con, sizeof(ULINT), 1, fp);
  size += sizeof(int) * fwrite(&MVPT->nodes_used, sizeof(int), 1, fp);
  size += sizeof(NODE) * fwrite(MVPT->nodes, sizeof(NODE), MVPT->nodes_used, fp);

  return size;  
}

int MVP_TREE_read(MVP_TREE *MVPT, FILE *fp)
{
  int size = 0;
  if (MVPT == NULL) return 0;

  if (MVPT->num_ref > 0) {
    MVPT->SPdist = mallocec(MVPT->num_ref * sizeof(DIST_TYPE));
  }

  size += sizeof(int) * fread(&MVPT->no_pts, sizeof(int), 1, fp);
  size += sizeof(ULINT) * fread(&MVPT->num_calc_con, sizeof(ULINT), 1, fp);
  size += sizeof(int) * fread(&MVPT->nodes_used, sizeof(int), 1, fp);

  MVPT->nodes_avail = MVPT->nodes_used;
  MVPT->nodes = mallocec(MVPT->nodes_used * sizeof(NODE));
  size += sizeof(NODE) * fread(MVPT->nodes, sizeof(NODE), MVPT->nodes_used, fp);

  return size;
}

char *MVP_TREE_sprint_stats(MVP_TREE *MVPT, int options)
{
  char *buf = NULL;
  int size = 0;
  int len = 0;
  absprintf(&buf, &size, &len, "\n*** MVP-tree Statistics ***\n");
  absprintf(&buf, &size, &len, "Data points: %d\n", MVPT->no_pts);
  absprintf(&buf, &size, &len, "Available nodes: %d\n", MVPT->nodes_avail);
  absprintf(&buf, &size, &len, "Used nodes: %d (%d bytes)\n", MVPT->nodes_used,
	    MVPT->nodes_used * sizeof(NODE));
  absprintf(&buf, &size, &len, "Distance evalutions for creation: %d\n", MVPT->num_calc_con);
  return buf;
}

void MVP_TREE_fprint_stats(MVP_TREE *MVPT, FILE *fp, int options) 
{
  char *buf = MVP_TREE_sprint_stats(MVPT, options);
  fprintf(fp, "%s", buf);
  free(buf);
}




void MVP_TREE_rng_srch(MVP_TREE *MVPT, BIOSEQ *query, int range)
{
  MVPT->seqs_visited = 0;
  MVPT->seqs_hit = 0;
  MVPT->nodes_visited = 1;
  MVPT->nodes_hit = 1;
  MVPT->hits_used = 0;
  search(MVPT, 1, 0, query, (DIST_TYPE) range);
}


/* ************************************************ */
/*  Fragment metric index - MVPT2                   */
/* ************************************************ */

static 
void get_pt(void *udb0, ULINT i, BIOSEQ *seq)
{
  UFRAG_DB *udb = (UFRAG_DB *) udb0;
  int offset = ufragdb_get_ufrag(udb, i);
  fastadb_get_Ffrag(udb->sdb, udb->frag_len, 
		    seq, offset);  
} 

static
int frag_metric(void *matrix, BIOSEQ *s1, BIOSEQ *s2)
{
  SCORE_MATRIX *D = (SCORE_MATRIX *) matrix;  
  return D->eval_score(D, s1->start, s2->start, s1->len);
}


MVP_TREE *MVP_TREE_init2(const char *db_name, ULINT frag_len, 
			 const char *matrix_name, int num_ref)
{
  UFRAG_DB *udb = ufragdb_init();
  SCORE_MATRIX *S;
  SCORE_MATRIX *D;
  MVP_TREE *MVPT2;
  DIST_TYPE **refdist = NULL;
  int i = 0;
  int *oid;
  int oid_size;

  ufragdb_create(udb, db_name, frag_len, 1);

  S = SCORE_MATRIX_from_file(matrix_name);
  S->set_conv_type(S, MAX);
  D = S->matrix_conv(S, NULL, 0);

  oid_size = ufragdb_get_nopts(udb);
  oid = mallocec(oid_size * sizeof(int));

  if (num_ref > 0) {
    refdist = mallocec(oid_size * sizeof(DIST_TYPE *));
    for (i=0; i < oid_size; i++)
      refdist[i] = callocec(num_ref, sizeof(DIST_TYPE));
  }

  MVPT2 = MVP_TREE_init(udb, D, get_pt, frag_metric, num_ref, refdist);

  for (i=0; i < oid_size; i++)
    oid[i] = i;

  printf("Before _create\n");
  MVP_TREE_create(MVPT2, oid, oid_size);
  MVP_TREE_fprint_stats(MVPT2, stdout, 0);

  SCORE_MATRIX_del(S);
  free(oid);
  
  return MVPT2;
} 

void MVP_TREE_del2(MVP_TREE *MVPT2)
{
  UFRAG_DB *udb = (UFRAG_DB *) MVPT2->db;
  SCORE_MATRIX *D = (SCORE_MATRIX *) MVPT2->matrix;
  int i;

  ufragdb_del(udb);
  SCORE_MATRIX_del(D);

  if (MVPT2->num_ref > 0) {
    for (i=0; i < MVPT2->no_pts; i++)
      free(MVPT2->refdist[i]);
    free(MVPT2->refdist);
  }

  MVP_TREE_destroy(MVPT2);
}


void MVP_TREE_save2(MVP_TREE *MVPT2, const char *filename) 
{
  FILE *fp = fopen(filename, "wb");
  UFRAG_DB *udb = (UFRAG_DB *) MVPT2->db;
  SCORE_MATRIX *D = (SCORE_MATRIX *) MVPT2->matrix;
  int i;

  if (fp == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "MVP_TREE_save2(): Could not open the file %s.",
		   filename);
  
  SCORE_MATRIX_write(D, fp);
  ufragdb_write(udb, fp);
  fwrite(&MVPT2->num_ref, sizeof(int), 1, fp);

  if (MVPT2->num_ref > 0) 
    for (i=0; i < MVPT2->no_pts; i++)
      fwrite(MVPT2->refdist[i], sizeof(DIST_TYPE), MVPT2->num_ref, fp);

  MVP_TREE_write(MVPT2, fp);
  fclose(fp);
}

MVP_TREE *MVP_TREE_load2(const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  UFRAG_DB *udb = ufragdb_init();
  SCORE_MATRIX *D;
  int i = 0;
  int num_ref = 0;
  DIST_TYPE **refdist = NULL;
  MVP_TREE *MVPT2;
  int oid_size;

  char *db_dir = NULL;
  char *db_base = NULL;


  if (fp == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "MVP_TREE_load2(): Could not open the file %s.",
		   filename);
  
  split_base_dir(filename, &db_base, &db_dir);

  D = SCORE_MATRIX_read(fp);
  ufragdb_read(udb, fp, db_dir);
  fread(&num_ref, sizeof(int), 1, fp);

  oid_size = ufragdb_get_nopts(udb);
  if (num_ref > 0) {
    refdist = mallocec(oid_size * sizeof(DIST_TYPE *));
    for (i=0; i < oid_size; i++) {
      refdist[i] = callocec(num_ref, sizeof(DIST_TYPE));
      fread(refdist[i], sizeof(DIST_TYPE), num_ref, fp);
    }
  }

  MVPT2 = MVP_TREE_init(udb, D, get_pt, frag_metric, num_ref, refdist);
  MVP_TREE_read(MVPT2, fp);

  fclose(fp);
  free(db_dir);
  free(db_base);

  return MVPT2;
}


HIT_LIST_t *MVP_TREE_rng_srch2(MVP_TREE *MVPT2, BIOSEQ *query, int range,
			       HIT_LIST_t *HL)
{
  UFRAG_DB *udb = (UFRAG_DB *) MVPT2->db;
  SCORE_MATRIX *D = (SCORE_MATRIX *) MVPT2->matrix;
  int j;
  BIOSEQ dupl;
  ULINT *dfrags;
  ULINT dsize;
  int k;
  int d0;

  if (HL == NULL) 
    HL = HIT_LIST_create(query, udb->sdb, D, range, -1, 0);
  else 
    HIT_LIST_reset(HL, query, udb->sdb, D, range, -1, 0);


  MVP_TREE_rng_srch(MVPT2, query, range);

  for (j=0; j < MVPT2->hits_used; j++) {
    get_pt(udb, MVPT2->hits[j], &dupl);
    d0 = (int) MVPT2->hit_dists[j];

    HIT_LIST_insert_hit(HL, &dupl, d0, 0); 
    HIT_LIST_count_seq_hit(HL, 1);
    
    ufragdb_get_dfrags2(udb, MVPT2->hits[j], &dfrags, &dsize);
    MVPT2->seqs_hit += dsize;
    for (k=0; k < dsize; k++) {
      fastadb_get_Ffrag(udb->sdb, udb->frag_len, &dupl, dfrags[k]);
      HIT_LIST_insert_hit(HL, &dupl, d0, 0); 
    }
  }
  HIT_LIST_count_FS_seq_visited(HL, MVPT2->nodes_visited);
  HIT_LIST_count_FS_seq_hit(HL, MVPT2->nodes_hit);
  HIT_LIST_count_seq_visited(HL, MVPT2->seqs_visited);


  return HL;
}


