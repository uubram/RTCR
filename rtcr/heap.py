#!/usr/bin/env python

class Heap:
    """Min heap implementation which allows for deleting specific items,
    identified by the item itself (using "==" not "is" operator) as opposed to
    its index in the heap, in logarithmic time.

    Note:
    - It uses a dictionary to keep track of the location of items in the heap.
    - Do not modify self.heap without also modifying self.inds.
    """
    def __init__(self):
        self.heap = []
        self.inds = {} # indexes of items in the heap

    def __len__(self):
        return len(self.heap) 

    def __contains__(self, item):
        return item in self.inds

    def _swap(self,i,j):
        """Swap index i and j."""
        self.heap[i], self.heap[j] = self.heap[j], self.heap[i]

        # set correct indices for the items
        self.inds[self.heap[i]] = i
        self.inds[self.heap[j]] = j

    def _bubble_up(self,i):
        parent = lambda i:max((i-1)//2,0)
        while self.heap[parent(i)] > self.heap[i]:
            self._swap(parent(i),i)
            i = parent(i)

    def _bubble_down(self,i):
        size = len(self.heap)
        while i < size:
            ci = 2*i + 1 # left child index
            ri = ci + 1 # right child index
            if (ci < size and self.heap[i] > self.heap[ci]) or \
                    (ri < size and \
                    self.heap[i] > self.heap[ri]):
                if ri < size and self.heap[ri] < self.heap[ci]:
                    ci = ri
                self._swap(i, ci) # swap with smallest child
                i = ci
            else:
                return

    def has_item(self, item):
        return item in self

    def has_heap_property(self):
        """Return True iff heap still has the heap property."""
        def check(i, size):
            li = 2*i + 1 # left child index
            ri = li + 1 # right child index
            if li < size:
                if self.heap[i] > self.heap[li]:
                    return False
                else:
                    return check(li, size)
            if ri < size:
                if self.heap[i] > self.heap[ri]:
                    return False
                else:
                    return check(ri, size)
            return True
        return check(0, len(self.heap))

    def peek(self):
        if len(self.heap) == 0:
            raise Exception("Cannot peek empty heap.")

        return self.heap[0]

    def push(self,item):
        """Push item onto the heap.

        :item: object that can be compared using "<" and ">" operators.
        """
        i = len(self.heap) # location of the new item
        if item in self:
            raise Exception("Item already on the heap.")
        self.heap += [item]
        self.inds[self.heap[i]] = i
        self._bubble_up(i) # restore the heap property

    def pop(self):
        """Pop minimum item from the heap."""
        if len(self.heap) == 0:
            raise Exception("Cannot pop from empty heap.")

        result = self.heap[0]
        self.remove_at(0)
        return result

    def remove_at(self,i):
        """Remove item at index i from the heap."""
        assert i < len(self.heap) # index out of range
        item = self.heap[i]

        # swap with right-most leaf and remove the item from
        # heap and inds.
        self._swap(i,-1)
        del self.heap[-1]
        del self.inds[item]

        # Handle edge cases
        if len(self.heap) == 0 or i >= len(self.heap):
            return

        # restore heap property
        parent = lambda i:max((i-1)//2,0)
        if self.heap[parent(i)] > self.heap[i]:
            self._bubble_up(i)
        else:
            self._bubble_down(i)

    def remove_item(self,item):
        """Remove specified item from the heap."""
        assert item in self.inds
        i = self.inds[item]
        assert item == self.heap[i] # are heap and inds still in sync
        self.remove_at(i)
