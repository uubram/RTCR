#include <Python.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include "util.h"

#ifndef PyMODINIT_FUNC /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

struct Pfm{
    uint32_t *freq;
    uint8_t *qual;
    ssize_t length;
};

typedef struct {
    PyObject_HEAD
    struct Pfm *pfm;            /* Position frequency matrix */
    char *seq;                  /* Consensus sequence */
    uint8_t *qual;              /* Quality score for every base in consensus */
    ssize_t length;             /* Length of seq (and of qual) */
    uint32_t count;             /* Number of sequences */
    uint32_t base_count;        /* Number of bases */
    uint32_t mutation_count;    /* Number of mutated bases */
} ConsensusQSequence;

#define MAX(A, B) ((A)>(B) ? (A):(B))
#define MIN(A, B) ((A)<(B) ? (A):(B))

static void
del_pfm(struct Pfm *a)
{
    if (a == NULL)
        return;

    if (a->freq != NULL)
        free(a->freq);
    if (a->qual != NULL)
        free(a->qual);
    free(a);
}

/*
 * Adds two position frequency matrices (pfms) together and returns the result
 * in a new pfm.
 *
 * @param offset: number of places (not positions!) pfm b is shifted to the
 * right of a. To calculate correct offset, multiply the number of positions b
 * is shifted to the right by the size of the alphabet. e.g.:
 * consensus of pfm a = "AACCTG"
 * consensus of pfm b =   "CCTG"
 * with alphabet = "ACGTN" correct offset would be 2 * 5 = 10.
 */
static struct Pfm *
pfm_add(struct Pfm *a, struct Pfm *b, ssize_t offset)
{
    if (offset < 0){
        struct Pfm *tmp = a;
        a = b;
        b = tmp;
        offset = -offset;
    }

    if (offset > a->length)
        return NULL;

    struct Pfm *c = safe_malloc(sizeof(*c));
    c->length = MAX(a->length, offset + b->length);
    c->freq = safe_malloc(sizeof(*c->freq) * c->length);
    c->qual = safe_malloc(sizeof(*c->qual) * c->length);

    int i;
    for (i = 0; i < c->length; i++){
        if (i < offset || i >= (offset + b->length)){
            c->freq[i] = a->freq[i];
            c->qual[i] = a->qual[i];
        }
        else if (i >= a->length){
            c->freq[i] = b->freq[i - offset];
            c->qual[i] = b->qual[i - offset];
        }
        else{
            c->freq[i] = a->freq[i] + b->freq[i - offset];
            c->qual[i] = MAX(a->qual[i], b->qual[i - offset]);
        }
    }
    return c;
}

static struct Pfm *
pfm_copy(const struct Pfm *a)
{
    struct Pfm *b = safe_malloc(sizeof(*b));
    b->length = a->length;
    b->freq = safe_malloc(sizeof(*b->freq) * b->length);
    b->qual = safe_malloc(sizeof(*b->qual) * b->length);
    memcpy(b->freq, a->freq, sizeof(*b->freq) * b->length);
    memcpy(b->qual, a->qual, sizeof(*b->qual) * b->length);
    return b;
}

/*
 * Updates several fields based on the pfm.
 *
 * Fields updated are: seq, qual, base_count, and mutation_count. Notably,
 * the 'count' field is *NOT* updated because that cannot be reliably inferred
 * from the pfm.
 */
static void
ConsensusQSequence_update(ConsensusQSequence *self){
    if (self->pfm == NULL){
        if (self->seq != NULL)
            free(self->seq);
        if (self->qual != NULL)
            free(self->qual);
        self->length = 0;
        self->count = 0;
        self->base_count = 0;
        self->mutation_count = 0;
        return;
    }

    const uint32_t length = self->pfm->length / 5;

    if (self->length != length || self->seq == NULL || self->qual == NULL){
        if (self->seq != NULL)
            free(self->seq);

        if (self->qual != NULL)
            free(self->qual);

        self->seq = safe_malloc(sizeof(*self->seq) * (length + 1));
        self->qual = safe_malloc(sizeof(*self->qual) * length);
    }

    self->length = length;
    self->base_count = 0;
    self->mutation_count = 0;

    ssize_t i;
    for (i = 0; i < self->length; i++){
        char cch = 'N'; /* consensus character */
        uint8_t cq = 0; /* consensus quality */
        uint32_t max_freq = 0;
        uint32_t base_count = 0;

        ssize_t j = 0;
        for (j = 0; j < 4; j++){ /* Get known base with highest freq */
            size_t pfm_pos = i*5 + j;
            base_count += self->pfm->freq[pfm_pos];
            if (self->pfm->freq[pfm_pos] > max_freq ||
                    (self->pfm->freq[pfm_pos] == max_freq &&
                     max_freq > 0 &&
                     self->pfm->qual[pfm_pos] > cq)){
                max_freq = self->pfm->freq[pfm_pos];
                cq = self->pfm->qual[pfm_pos];
                switch(j){
                case 0:
                    cch = 'A';
                    break;
                case 1:
                    cch = 'C';
                    break;
                case 2:
                    cch = 'G';
                    break;
                case 3:
                    cch = 'T';
                    break;
                }
            }
            else if (self->pfm->freq[pfm_pos] == max_freq && max_freq > 0 &&
                    self->pfm->qual[pfm_pos] == cq){
                cch = 'N';
                /* Note, do not put cq to 0 here as it is still needed to find
                 * out if another base with the same freq has a higher quality
                 * score. */
            }
        }
        base_count += self->pfm->freq[i * 5 + 4];
        if (cch == 'N')
            cq = 0;
        if (max_freq == 0){
            max_freq = self->pfm->freq[i * 5 + 4];
            cch = 'N';
            cq = 0;
        }
        self->seq[i] = cch;
        self->qual[i] = cq;
        self->base_count += base_count;
        self->mutation_count += base_count - max_freq;
    }
    self->seq[self->length] = '\0';
}

static void
ConsensusQSequence_dealloc(ConsensusQSequence *self)
{
    if (self->seq != NULL)
        free(self->seq);

    if (self->qual != NULL)
        free(self->qual);

    if (self->pfm != NULL)
        del_pfm(self->pfm);

    self->ob_type->tp_free((PyObject *)self);
}

static PyObject *
ConsensusQSequence_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    (void)args;
    (void)kwds;
    ConsensusQSequence *self;

    self = (ConsensusQSequence *)type->tp_alloc(type, 0);
    if (self != NULL){
        self->pfm = NULL;
        self->seq = NULL;
        self->qual = NULL;
        self->length = 0;
        self->count = 0;
        self->base_count = 0;
        self->mutation_count = 0;
    }

    return (PyObject *)self;
}

static PyTypeObject ConsensusQSequenceType;

static int
ConsensusQSequence_init(ConsensusQSequence *self, PyObject *args,
        PyObject *kwds)
{
    PyObject *other = NULL;
    self->count = 0;

    static char *kwlist[] = {"other", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &other))
        return -1;

    if (other == NULL || other == Py_None ||
            (PySequence_Check(other) != 0 && PySequence_Length(other) == 0)){
        return 0;
    }

    if (PyObject_TypeCheck(other, &ConsensusQSequenceType) != 0){
        ConsensusQSequence *other_cs = (ConsensusQSequence *)other;
        if (other_cs->pfm != NULL){
            self->pfm = pfm_copy(other_cs->pfm);
            self->count = other_cs->count;
            /*ConsensusQSequence_update(self);*/
            const size_t length = self->pfm->length / 5;
            self->seq = safe_malloc(sizeof(*self->seq) * (length + 1));
            self->qual = safe_malloc(sizeof(*self->qual) * length);
            memcpy(self->seq, other_cs->seq,
                    sizeof(*self->seq) * (length+1));
            memcpy(self->qual, other_cs->qual,
                    sizeof(*self->qual) * length);
            self->length = other_cs->length;
            self->base_count = other_cs->base_count;
            self->mutation_count = other_cs->mutation_count;
        }
        return 0;
    }

    /* In an unpickling scenario, the argument will be a tuple containing
     * two other tuples (freq and qual fields of a pfm respectively) and a
     * count representing the number of sequences. */
    if (PyTuple_Check(other) != 0 && PyTuple_GET_SIZE(other) == 3 &&
            PyTuple_Check(PyTuple_GET_ITEM(other, 0))){
        PyObject *freq = PyTuple_GET_ITEM(other, 0);
        PyObject *qual = PyTuple_GET_ITEM(other, 1);
        PyObject *count = PyTuple_GET_ITEM(other, 2);
        
        self->count = PyInt_AsLong(count);
        if (PyErr_Occurred() != NULL)
            return -1; 

        Py_ssize_t size = PyTuple_GET_SIZE(freq);
        if (size != PyTuple_GET_SIZE(qual)){
            PyErr_SetString(PyExc_ValueError,
                    "freq and qual in tuple (first two items) are not the "
                    "same length");
            return -1;
        }

        if (size % 5 != 0){
            PyErr_SetString(PyExc_ValueError,
                    "length of freq (first item in tuple) is not "
                    "a multiple of 5.");
            return -1;
        }

        self->pfm = safe_malloc(sizeof(*self->pfm));
        self->pfm->length = size;
        self->pfm->freq = safe_calloc(self->pfm->length,
                sizeof(*self->pfm->freq));
        self->pfm->qual = safe_calloc(self->pfm->length,
                sizeof(*self->pfm->qual));
        for (ssize_t i = 0; i < self->pfm->length; i++){
            self->pfm->freq[i] = PyInt_AsLong(PyTuple_GET_ITEM(freq, i));
            self->pfm->qual[i] = PyInt_AsLong(PyTuple_GET_ITEM(qual, i));
            if (PyErr_Occurred() != NULL)
                return -1;
        }
        ConsensusQSequence_update(self);
        if (self->base_count < self->count){
            PyErr_SetString(PyExc_ValueError,
                    "number of bases in pfm exceeds sequence count.");
            return -1;
        }
        return 0;
    }

    PyObject *seq = NULL;
    PyObject *qual = NULL;
    PyObject *count = NULL;
    PyObject *item = NULL;

    if (PySequence_Check(other) != 0){
        seq = PySequence_GetItem(other, 0);
        qual = PySequence_GetItem(other, 1);
        if (PyErr_Occurred() != NULL)
            goto cleanup;
        count = PySequence_GetItem(other, 2);
        PyErr_Clear();
    }
    else{
        seq = PyObject_GetAttrString(other, "seq");
        qual = PyObject_GetAttrString(other, "qual");
        if (PyErr_Occurred() != NULL)
            goto cleanup;
        count = PyObject_GetAttrString(other, "count");
        PyErr_Clear();
    }

    if (PyString_Check(seq) == 0){
        PyErr_SetString(PyExc_TypeError, "seq is not a string.");
        goto cleanup;
    }

    if (PySequence_Check(qual) == 0){
        PyErr_SetString(PyExc_TypeError,
                "qual does not provide sequence protocol "
                "(i.e. it is not a list or tuple-like object).");
        goto cleanup;
    }

    if (count != NULL){
        self->count = PyInt_AsLong(count);
        if (PyErr_Occurred() != NULL)
            goto cleanup;
    }
    else {
        self->count = 1;
    }

    if (self->count < 1){
        PyErr_SetString(PyExc_ValueError, "count < 1");
        goto cleanup;
    }

    char *s = PyString_AsString(seq);
    if (PyErr_Occurred() != NULL)
        goto cleanup;
    
    self->length = strlen(s);

    if (self->length < 1){
        PyErr_SetString(PyExc_ValueError, "len(seq) < 1");
        goto cleanup;
    }

    if (PySequence_Size(qual) != self->length){
        PyErr_SetString(PyExc_ValueError, "len(seq) != len(qual)");
        goto cleanup;
    }

    self->pfm = safe_malloc(sizeof(*self->pfm));
    self->pfm->length = self->length * 5;
    self->pfm->freq = safe_calloc(self->pfm->length, sizeof(*self->pfm->freq));
    self->pfm->qual = safe_calloc(self->pfm->length, sizeof(*self->pfm->qual));

    Py_ssize_t i;
    for (i = 0; i < self->length; i++){
        item = PySequence_GetItem(qual, i);
        if (item == NULL)
            goto cleanup;
        long Q = PyInt_AsLong(item);
        if (PyErr_Occurred() != NULL)
            goto cleanup;
        int j;
        switch(*(s + i)){
        case 'A':
            j = 0;
            break;
        case 'C':
            j = 1;
            break;
        case 'G':
            j = 2;
            break;
        case 'T':
            j = 3;
            break;
        case 'N':
            j = 4;
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "Wrong character in seq.");
            goto cleanup;
        }

        int pfm_pos = i * 5 + j;
        self->pfm->freq[pfm_pos] = self->count;
        if (Q < 0 || Q > 41){
            PyErr_SetString(PyExc_ValueError,
                    "Quality not a valid PHRED score.");
            goto cleanup;
        }
        if (Q != 0 && j == 4){
            PyErr_SetString(PyExc_ValueError,
                    "Quality score of \'N\' should be 0.");
            goto cleanup;
        }
        self->pfm->qual[pfm_pos] = Q;
    }

    ConsensusQSequence_update(self);
    Py_XDECREF(seq);
    Py_XDECREF(qual);
    Py_XDECREF(count);
    Py_XDECREF(item);
    return 0;

    cleanup:
        Py_XDECREF(seq);
        Py_XDECREF(qual);
        Py_XDECREF(count);
        Py_XDECREF(item);
        return -1;
}

static PyObject *
ConsensusQSequence_add(ConsensusQSequence *self, PyObject *args,
        PyObject *kwds)
{
    int offset = 0;
    ConsensusQSequence *other;

    static char *kwlist[] = {"other", "offset", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist,
                &ConsensusQSequenceType, &other, &offset))
        return NULL;

    if (self->pfm != NULL && other->pfm != NULL){
        struct Pfm* c = pfm_add(self->pfm, other->pfm, offset * 5);
        if (c == NULL){
            PyErr_SetString(PyExc_ValueError,
                    "Unable to add position frequency matrices "
                    "(offset wrong?)");
            return NULL;
        }
        del_pfm(self->pfm);
        self->pfm = c;
        self->count += other->count;
    }
    else if (self->pfm == NULL){
        self->pfm = pfm_copy(other->pfm);
        self->count = other->count;
    }
    else{
        Py_RETURN_NONE;
    }

    ConsensusQSequence_update(self);
    Py_RETURN_NONE;
}

static PyObject *
ConsensusQSequence_reduce(ConsensusQSequence *self)
{
    if (self->pfm != NULL){
        PyObject *freq = PyTuple_New(self->pfm->length);
        PyObject *qual = PyTuple_New(self->pfm->length);
        PyObject *count = PyInt_FromLong(self->count);

        for (int i = 0; i < self->pfm->length; i++){
            PyTuple_SET_ITEM(freq, i, PyInt_FromLong(self->pfm->freq[i]));
            PyTuple_SET_ITEM(qual, i, PyInt_FromLong(self->pfm->qual[i]));
        }
        return PyTuple_Pack(2, Py_TYPE(self),
                PyTuple_Pack(1, PyTuple_Pack(3, freq, qual, count)));
    }
    else{
        return PyTuple_Pack(2, Py_TYPE(self), PyTuple_New(0));
    }
}

static PyObject *
ConsensusQSequence_get_seq(ConsensusQSequence *self, void *closure)
{
    (void)closure;
    if (self->seq == NULL)
        return Py_BuildValue("s", &"");
    else
        return Py_BuildValue("s", self->seq);
}

static PyObject *
ConsensusQSequence_get_qual(ConsensusQSequence *self, void *closure)
{
    (void)closure;
    PyObject *qual = PyList_New(self->length);
    Py_ssize_t i;
    for (i = 0; i < self->length; i++)
        PyList_SET_ITEM(qual, i, Py_BuildValue("i", self->qual[i]));
    return qual;
}

static PyObject *
ConsensusQSequence_get_count(ConsensusQSequence *self, void *closure)
{
    (void)closure;
    return Py_BuildValue("i", self->count);
}

static PyObject *
ConsensusQSequence_get_base_count(ConsensusQSequence *self, void *closure)
{
    (void)closure;
    return Py_BuildValue("i", self->base_count);
}

static PyObject *
ConsensusQSequence_get_mutation_count(ConsensusQSequence *self, void *closure)
{
    (void)closure;
    return Py_BuildValue("i", self->mutation_count);
}

static PyMethodDef ConsensusQSequence_methods[] = {
    {"add", (PyCFunction)ConsensusQSequence_add,
        METH_VARARGS | METH_KEYWORDS, ""},
    {"__reduce__", (PyCFunction)ConsensusQSequence_reduce,
        METH_NOARGS, ""},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static PyGetSetDef ConsensusQSequence_getset[] = {
    {"seq", (getter)ConsensusQSequence_get_seq,
        NULL, NULL, NULL},
    {"qual", (getter)ConsensusQSequence_get_qual,
        NULL, NULL, NULL},
    {"count", (getter)ConsensusQSequence_get_count,
        NULL, NULL, NULL},
    {"base_count", (getter)ConsensusQSequence_get_base_count,
        NULL, NULL, NULL},
    {"mutation_count", (getter)ConsensusQSequence_get_mutation_count,
        NULL, NULL, NULL},
    {NULL, NULL, NULL, NULL, NULL} /* Sentinel */
};

static PyTypeObject ConsensusQSequenceType = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /* ob_size */
    "cseq.ConsensusQSequence",             /* tp_name */
    sizeof(ConsensusQSequence),                 /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)ConsensusQSequence_dealloc,     /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_compare */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    0,                                          /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    ConsensusQSequence_methods,                 /* tp_methods */
    0,                                          /* tp_members */
    ConsensusQSequence_getset,                  /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)ConsensusQSequence_init,          /* tp_init */
    0,                                          /* tp_alloc */
    ConsensusQSequence_new,                     /* tp_new */
    0,                                          /* tp_free */
    0,                                          /* tp_is_gc */
    0,                                          /* tp_bases */
    0,                                          /* tp_mro */
    0,                                          /* tp_cache */
    0,                                          /* tp_subclasses */
    0,                                          /* tp_weaklist */
    0,                                          /* tp_del */
    0                                           /* tp_version_tag */
};

static PyObject *
merge_stats(PyObject *self, PyObject *args, PyObject *kwds)
{
    (void)self;
    ConsensusQSequence *a = NULL;
    ConsensusQSequence *b = NULL;
    int offset = 0;

    static char *kwlist[] = {"a", "b", "offset", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|i", kwlist,
                &ConsensusQSequenceType, &a,
                &ConsensusQSequenceType, &b,
                &offset))
        return NULL;
   
    if (offset < 0){
        ConsensusQSequence *tmp = a;
        a = b;
        b = tmp;
        offset = -offset;
    }

    long hd = 0;
    long merge_score = 0;
    long max_mut_q = 0;

    if (a->length > 0 && b->length > 0){
        int max_i = MIN(a->length, b->length + offset);
        for (int i = offset; i < max_i; i++){
            if (a->seq[i] != b->seq[i - offset]){
                int minq = MIN(a->qual[i], b->qual[i - offset]); 
                merge_score += minq;
                max_mut_q = MAX(max_mut_q, minq);
                hd++;
            }
        }
    }

    return PyTuple_Pack(3, PyInt_FromLong(hd), PyInt_FromLong(merge_score),
                PyInt_FromLong(max_mut_q));
}

static int
score_matches(const char *s1, const char *s2, int len)
{
    int score = 0;
    for (int i = 0; i < len; i++){
        if (s1[i] == 'N' || s2[i] == 'N')
            continue;
        if (s1[i] == s2[i])
            score += 1;
        else
            score -= 1;
    }
    return score;
}

static PyObject *
get_offset(PyObject *self, PyObject *args, PyObject *kwds)
{
    (void)self;
    char *s1 = NULL;
    char *s2 = NULL;
    int k;
    int moff;

    static char *kwlist[] = {"s1", "s2", "k", "moff", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ssii", kwlist, &s1, &s2, &k,
                &moff))
        return NULL;

    int len_s1 = strlen(s1);
    int len_s2 = strlen(s2);

    /* Make sure s2 is the longer of the two strings because different kmers
     * will be selected within s2 */
    int swapped = 0;
    if (len_s1 > len_s2){
        char *tmp_s = s1;
        int tmp_len = len_s1;
        s1 = s2;
        len_s1 = len_s2;
        s2 = tmp_s;
        len_s2 = tmp_len;
        swapped = 1;
    }

    int mid = len_s1 / 2;

    if (k > MIN(len_s1, len_s2)){
        PyErr_SetString(PyExc_ValueError, "String(s) shorter than window");
        return NULL;
    }

    int start = mid - (k / 2);
    int moff_left = MIN(moff, start);
    int moff_right = MIN(moff, len_s2 - (start + k));

    const int cur_moff = MAX(moff_left, moff_right);
    const int NO_OFFSET = moff + 1;

    int best_offset = NO_OFFSET;
    int best_m = 0;
    for (int offset = 0; offset <= cur_moff; offset++){
        int m = 0;
        if (offset <= moff_right){
            m = score_matches(s1 + start, s2 + start + offset, k);
            if (m == k)
                return swapped?
                    PyTuple_Pack(2, PyInt_FromLong(-offset), PyInt_FromLong(m))
                    :
                    PyTuple_Pack(2, PyInt_FromLong(offset), PyInt_FromLong(m));
        }
        if (best_offset == NO_OFFSET || best_m < m){
            best_offset = offset;
            best_m = m;
        }
        if (offset <= moff_left){
            m = score_matches(s1 + start, s2 + start - offset, k);
            if (m == k)
                return swapped?
                    PyTuple_Pack(2, PyInt_FromLong(offset), PyInt_FromLong(m))
                    :
                    PyTuple_Pack(2, PyInt_FromLong(-offset),
                        PyInt_FromLong(m));
        }
        if (best_offset == NO_OFFSET || best_m < m){
            best_offset = -offset;
            best_m = m;
        }
    }

    if (best_offset == NO_OFFSET)
        Py_RETURN_NONE;
    else
        return swapped?
            PyTuple_Pack(2, PyInt_FromLong(-best_offset),
                    PyInt_FromLong(best_m))
            :
            PyTuple_Pack(2, PyInt_FromLong(best_offset),
                    PyInt_FromLong(best_m));
}

/*
 * Find start and end of the alignment intervals on the query and
 * reference sequences, using the CIGAR string. Returns 4-tuple
 * (qas, qae, ras, rae), where qas and qae are the 0-based alignment
 * start- and end-positions on the query respectively (SEQ field of a SAM
 * record), followed by a similar pair for the reference.
 */
static PyObject *
cigar_intervals(PyObject *self, PyObject *args, PyObject *kwds)
{
    (void)self;

    const char *cigar;
    int ras; /* 0-based(!) leftmost mapping position on reference.*/
    /* Note, the 'POS' field of a SAM record is 1-based, therefore subtract 1
     * from that field when passing it to this method.*/

    static char *kwlist[] = {"cigar", "ras", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "si", kwlist, (char *)&cigar,
                &ras))
        return NULL;

    int rae = ras;
    int qas = 0;
    int qae = qas;

    int charnumber = 0; /* number of non-numeric cigar characters processed */
    bool hardclipped_start = false;
    while (*cigar != '\0'){
        int n = 0;
        while(isdigit(*cigar)){
             n = n*10 + (*cigar - '0');
             cigar++;
        }
        switch(*cigar){
        case 'M':
        case 'X':
        case '=':
            rae += n;
            qae += n;
            break;
        case 'I': /* insertion to the reference */
            qae += n;
            break;
        case 'D': /* deletion from the reference */
            rae += n;
            break;
        case 'N': /* skipped region from the reference */
            rae += n;
            break;
        case 'S': /* soft clip on the query */
            if (charnumber == 0 || (charnumber == 1 && hardclipped_start)){
                qas += n;
                qae += n;
            }
            break;
        case 'H':
            if (charnumber == 0)
                hardclipped_start = true;
            break;
        case 'P':
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "non-cigar character");
            return NULL;
        }
        if (*cigar != '\0')
            cigar++;
        charnumber++;
    }
    return PyTuple_Pack(4, PyInt_FromLong(qas), PyInt_FromLong(qae),
           PyInt_FromLong(ras), PyInt_FromLong(rae));
}

/*Find which position on the query corresponds to the given position
 * on the reference.
 */
static PyObject *
cigar_rpos2qpos(PyObject *self, PyObject *args, PyObject *kwds)
{
    (void)self;

    const char *cigar; /* cigar string from a SAM record */
    int ras; /* 0-based(!) leftmost mapping position on reference.*/
    /* Note, the 'POS' field of a SAM record is 1-based, therefore subtract 1
     * from that field when passing it to this method.*/
    int target_rpos; /* target position on the reference sequence */

    static char *kwlist[] = {"cigar", "ras", "target_rpos", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sii", kwlist, (char *)&cigar,
                &ras, &target_rpos))
        return NULL;

    int rpos = ras;
    int qas = 0;
    int qpos = qas;
    
    if (target_rpos < ras)
        Py_RETURN_NONE;

    while (*cigar != '\0'){
        int n = 0;
        while(isdigit(*cigar)){
             n = n*10 + (*cigar - '0');
             cigar++;
        }
        switch(*cigar){
        case 'M':
        case 'X':
        case '=':
            rpos += n;
            qpos += n;
            if (target_rpos < rpos)
                return PyInt_FromLong(qpos - (rpos - target_rpos));
            break;
        case 'I': /* insertion to the reference */
            qpos += n;
            break;
        case 'D': /* deletion from the reference */
            rpos += n;
            break;
        case 'N': /* skipped region from the reference */
            rpos += n;
            break;
        case 'S': /* soft clip on the query */
            qpos += n;
            break;
        case 'H':
        case 'P':
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "non-cigar character");
            return NULL;
        }
        if (*cigar != '\0')
            cigar++;
        if (target_rpos < rpos)
            Py_RETURN_NONE;
    }
    Py_RETURN_NONE;
}

static PyObject *
get_error_stats(PyObject *self, PyObject *args, PyObject *kwds)
{
    (void)self;

    PyObject *aln;
    char *rseq;
    PyObject *Q_mm;
    PyObject *Q_n;
    int r_roi_start = 0;
    int r_roi_end = 0;

    static char *kwlist[] = {"aln", "rseq", "Q_mm", "Q_n",
        "r_roi_start", "r_roi_end", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OsO!O!|ii", kwlist, &aln,
                &rseq,
                &PyList_Type, &Q_mm,
                &PyList_Type, &Q_n,
                &r_roi_start, &r_roi_end))
        return NULL;

    if (r_roi_end == 0)
        r_roi_end = strlen(rseq);

    if (r_roi_start < 0){
        PyErr_SetString(PyExc_ValueError, "roi starting before position 0");
        return NULL;
    }

    if (r_roi_start >= r_roi_end){
        PyErr_SetString(PyExc_ValueError, "roi starts after its end");
        return NULL;
    }

    PyObject *qseq_pyobj = PyObject_GetAttrString(aln, "SEQ");
    PyObject *qualstr_pyobj = PyObject_GetAttrString(aln, "QUAL");
    PyObject *cigar_pyobj = PyObject_GetAttrString(aln, "CIGAR");
    PyObject *rpos_pyobj = PyObject_GetAttrString(aln, "POS");

    if (PyErr_Occurred() != NULL)
        goto cleanup;

    const char *qseq = PyString_AsString(qseq_pyobj);
    if (qseq == NULL)
        goto cleanup;

    /* ascii encoded (Phred + 33) quality string */
    const char *qualstr = PyString_AsString(qualstr_pyobj);
    if (qualstr == NULL)
        goto cleanup;

    const char *cigar = PyString_AsString(cigar_pyobj);
    if (cigar == NULL)
        goto cleanup;

    long rpos = PyInt_AsLong(rpos_pyobj) - 1;
    if (PyErr_Occurred() != NULL)
        goto cleanup;

    int qlen = 0; /* length of query within aligned part of the ROI */
    int mm = 0;
    int ins = 0;
    int dels = 0;

    int qpos = 0; /* position on query */

    int r_roi_as; /* alignment start index within ROI */
    if (rpos > r_roi_start)
        r_roi_as = rpos - r_roi_start;
    else
        r_roi_as = 0;
    int r_roi_ae = r_roi_as; /* alignment end index within ROI */

    while (*cigar != '\0'){
        if (rpos >= r_roi_end)
            break;
        int n = 0;
        while(isdigit(*cigar)){
             n = n*10 + (*cigar - '0');
             cigar++;
        }
        switch(*cigar){
        case 'M':
        case 'X':
        case '=':
            for (int j = 0; j < n; j++, rpos++, qpos++){
                if (rpos < r_roi_start)
                    continue;
                if (rpos >= r_roi_end)
                    break;
                r_roi_ae++;
                qlen++;
                int Q = *(qualstr + qpos) - 33;
                PyObject *elem = PyList_GetItem(Q_n, Q);
                if (elem == NULL)
                    goto cleanup;
                long val = PyInt_AsLong(elem);
                if (PyErr_Occurred() != NULL)
                    goto cleanup;
                PyList_SetItem(Q_n, Q, PyInt_FromLong(val + 1));
                if (*(rseq + rpos) != *(qseq + qpos)){
                    mm++;
                    PyObject *elem = PyList_GetItem(Q_mm, Q);
                    if (elem == NULL)
                        goto cleanup;
                    long val = PyInt_AsLong(elem);
                    if (PyErr_Occurred() != NULL)
                        goto cleanup;
                    PyList_SetItem(Q_mm, Q, PyInt_FromLong(val + 1));
                }
            }
            break;
        case 'I': /* insertion to the reference */
            qpos += n;
            if (rpos > r_roi_start){
                ins += n;
                qlen += n;
                for (int j = 0; j < n; j++){
                    int Q = *(qualstr + qpos) - 33;
                    PyObject *elem = PyList_GetItem(Q_n, Q);
                    if (elem == NULL)
                        goto cleanup;
                    long val = PyInt_AsLong(elem);
                    if (PyErr_Occurred() != NULL)
                        goto cleanup;
                    PyList_SetItem(Q_n, Q, PyInt_FromLong(val + 1));
                }
            }
            break;
        case 'D': /* deletion from the reference */
            for (int j = 0; j < n; j++, rpos++){
                if (rpos < r_roi_start)
                    continue;

                if (rpos > r_roi_start && rpos < r_roi_end){
                    dels++;
                    r_roi_ae++;
                }
            }
            break;
        case 'N': /* skipped region from the reference */
            rpos += n;
            break;
        case 'S': /* soft clip on the query */
            qpos += n;
            break;
        case 'H':
        case 'P':
            break;
        default:
            PyErr_SetString(PyExc_ValueError, "non-cigar character");
            goto cleanup;
        }
        if (*cigar != '\0')
            cigar++;
    }

    Py_XDECREF(qseq_pyobj);
    Py_XDECREF(qualstr_pyobj);
    Py_XDECREF(cigar_pyobj);
    Py_XDECREF(rpos_pyobj);

    return PyTuple_Pack(6, PyInt_FromLong(qlen), PyInt_FromLong(mm),
           PyInt_FromLong(ins), PyInt_FromLong(dels),
           PyInt_FromLong(r_roi_as), PyInt_FromLong(r_roi_ae));

    cleanup:
        Py_XDECREF(qseq_pyobj);
        Py_XDECREF(qualstr_pyobj);
        Py_XDECREF(cigar_pyobj);
        Py_XDECREF(rpos_pyobj);
        return NULL;
}

static PyMethodDef methods[] = {
    {"merge_stats", (PyCFunction)merge_stats,
        METH_VARARGS | METH_KEYWORDS, NULL},
    {"get_offset", (PyCFunction)get_offset,
        METH_VARARGS | METH_KEYWORDS, NULL},
    {"cigar_intervals", (PyCFunction)cigar_intervals,
        METH_VARARGS | METH_KEYWORDS, NULL},
    {"cigar_rpos2qpos", (PyCFunction)cigar_rpos2qpos,
        METH_VARARGS | METH_KEYWORDS, NULL},
    {"get_error_stats", (PyCFunction)get_error_stats,
        METH_VARARGS | METH_KEYWORDS, NULL},
    {NULL, NULL, 0, NULL}, /* Sentinel */
};

PyMODINIT_FUNC
initcseq(void)
{
    PyObject* m;

    if (PyType_Ready(&ConsensusQSequenceType) < 0)
        return;

    m = Py_InitModule("cseq", methods);

    if (m == NULL)
        return;

    /* Using precasted ConsensusQSequenceType prevents type-punning error when
     * compiling with strict-aliasing.
     */
    PyObject *op = (PyObject*) &ConsensusQSequenceType;
    Py_INCREF(op);

    PyModule_AddObject(m, "ConsensusQSequence",
            (PyObject *)&ConsensusQSequenceType);
}
