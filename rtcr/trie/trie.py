class Node(object):
    def __init__(self, parent, char, value):
        self._parent = parent
        self._char = char
        self._children = {}
        self._value = value
        self._isset = False

    def add_child(self, ch, value):
        if ch in self._children:
            return self._children[ch]

        node = Node(parent = self, char = ch, value = value)
        self._children[node.char] = node
        return node

    def remove_child(self, char):
        if char in self._children:
            del self._children[char]

    def is_leaf(self):
        return len(self._children) == 0

    def has_children(self):
        return len(self._children) > 0

    @property
    def parent(self):
        return self._parent

    @property
    def char(self):
        return self._char

    @property
    def children(self):
        return self._children

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        self._isset = True
        self._value = value

    @property
    def isset(self):
        return self._isset

    def clear(self):
        self._value = None
        self._isset = False

    def __repr__(self):
        return "%s(parent: %s, char: \"%s\", children: %s, value: %r)"%(
                self.__class__,
                str(self._parent),
                self._char,
                self._children.keys(),
                self._value)
    
    def __str__(self):
        return "%s -> %s"%(self._char, self._children.keys())

class Trie(object):
    """Trie structure, allowing strings to be inserted and associated with a
    value.
    """

    def __init__(self):
        self._root = Node(parent = None, char = None, value = None)
        self._length = 0

    def get_node(self, s):
        """Get node corresponding to the given string.
        
        None is returned if node does not exist.
        """
        node = self._root
        for ch in s:
            if not ch in node.children:
                return None
            else:
                node = node.children[ch]
        return node

    def _add(self, s):
        """Add string to the trie."""
        if not isinstance(s, basestring):
            raise TypeError("String is not (derived from) basestring.")

        node = self._root
        for ch in s:
            node = node.add_child(ch, value = None)
        return node

    def has_node(self, key):
        """Return True if a prefix exists in the trie.
        
        Note, this will also return true for nodes not associated with a
        particular value (i.e. a string corresponding to prefix has not
        explicitly been inserted into the trie).
        """
        node = self._root
        for ch in key:
            if not ch in node.children:
                return False
            else:
                node = node.children[ch]
        return True

    def __len__(self):
        return self._length

    def setdefault(self, key, defval = None):
        if key in self:
            return self[key]
        else:
            self[key] = defval
            return defval

    def __setitem__(self, key, value):
        """Insert a string into the trie and associate it with the given value
        (can be anything, including None).
        
        A TypeError will be raised if the string is not derived from
        basestring.
        """

        node = self._add(key)
        if not node.isset:
            self._length += 1
        node.value = value

    def __getitem__(self, key):
        """Get value associated with the given string from the trie.
        
        Note, None is considered a value. An exception will be thrown if string
        is not associated with a value.
        """
        if not isinstance(key, basestring):
            raise TypeError("key is not a string")
        node = self._root
        len_key = len(key)
        if len_key == 0 and self._root.isset:
            return self._root.value

        i = 0
        for i, ch in enumerate(key):
            if not ch in node.children:
                raise KeyError("String does not exist.")
            else:
                node = node.children[ch]
        if i+1 != len_key or not node.isset:
            raise KeyError("String does not exist.")
        return node.value

    def get(self, key, default = None):
        try:
            return self[key]
        except:
            return default

    def __contains__(self, key):
        try:
            self[key]
            return True
        except KeyError:
            return False

    def has_key(self, key):
        return key in self

    def __delitem__(self, key):
        """Remove all nodes that lead up to and only to the given string.

        If the string itself is a prefix of another string, then the string
        will be given a value of None (and no nodes will be removed).
        """
        node = self.get_node(key)

        if node is None or not node.isset:
            raise ValueError("String (\"%s\") does not exist."%s)

        node.clear()

        while node.is_leaf() and not node.isset and not node.parent is None:
            node.parent.remove_child(node.char)
            node = node.parent
        self._length -= 1
    
    def get_nearest_variants(self, s, maxhd = None):
        """Find string or strings most similar to s, but not s itself.
        
        Returns generator, yielding tuples of Hamming distance, string,
        and value corresponding to found string.

        :maxhd: maximum Hamming distance the function will consider for finding
        nearest variant to s.
        """
        stack = [(self._root, "", 0, 0)]
        next_stack = []
        cur_maxhd = 1
        done = False
        while stack:
            node, alt_s, pos, hd = stack.pop()

            if pos == len(s) and hd > 0 and node.isset:
                done = True
                yield hd, alt_s, node.value
                continue

            for ch, next_node in node.children.iteritems():
                if ch == s[pos]:
                    stack.append((next_node, alt_s + ch, pos + 1, hd))
                elif hd < cur_maxhd:
                    stack.append((next_node, alt_s + ch, pos + 1, hd + 1))
                elif not done:
                    # This next_node is too far, do not visit it unless nothing
                    # is found at current maxhd.
                    next_stack.append((next_node, alt_s + ch, pos + 1, hd + 1))

            if not stack and not done:
                if cur_maxhd == maxhd:
                    return
                cur_maxhd += 1
                stack = next_stack
                next_stack = []

    def neighbors(self, s, maxhd):
        """Search for all strings within the given Hamming distance of the
        given string.
        
        Returns generator, yielding tuples of Hamming distance, string,
        and value corresponding to found string.
        """
        if maxhd < 1:
            raise ValueError("maxhd < 1")
        if self.get_node(s) is None:
            raise ValueError("String (%s) does not exist."%s)
        stack = [(self._root, "", 0, 0)]
        while stack:
            node, alt_s, pos, hd = stack.pop()

            if pos == len(s) and node.isset and hd > 0:
                yield hd, alt_s, node.value
                continue

            for ch, next_node in node.children.iteritems():
                if ch == s[pos]:
                    stack.append((next_node, alt_s + ch, pos + 1, hd))
                elif hd < maxhd:
                    stack.append((next_node, alt_s + ch, pos + 1, hd + 1))

    def pairs(self, keylen, maxhd):
        """Generator function to iterate all pairs of strings of given length
        that are within a certain Hamming distance of each other.

        Yields tuples with the following format:
        (Hamming distance, string1, value1, string2, value2)
        """
        if keylen < 1:
            raise ValueError("keylen < 1")
        if maxhd < 1:
            raise ValueError("maxhd < 1")
        # Traverse the tree in search of strings of length keylen.
        targets = {} 
        stack = [(self._root, "", 0)]
        while stack:
            node, s, depth = stack.pop()

            if depth == keylen and node.isset:
                targets[s] = node

            for ch, next_node in node.children.iteritems():
                stack.append((next_node, s + ch, depth + 1))

        explored = set()
        for s, node_s in targets.iteritems():
            explored.add(s)
            stack = [(self._root, "", 0, 0)]
            while stack:
                node, alt_s, pos, hd = stack.pop()

                if pos == keylen:
                    if alt_s in targets:
                        yield hd, s, node_s.value, \
                                alt_s, targets[alt_s].value
                    # No need to go deeper into the tree because only equal
                    # length strings are considered.
                    continue

                n_explored = 0
                for ch, next_node in node.children.iteritems():
                    next_alt_s = alt_s + ch
                    if next_alt_s in explored:
                        # all pairs within maxhd for len(s) prefixes have been
                        # found already for the next_alt_s branch in the tree.
                        n_explored += 1
                        continue
                    if ch == s[pos]:
                        stack.append((next_node, next_alt_s, pos + 1, hd))
                    elif hd < maxhd:
                        stack.append((next_node, next_alt_s, pos + 1, hd + 1))

                if n_explored == len(node.children):
                    # Closing off current branch of the tree for current
                    # sequence length.
                    explored.add(alt_s)
                    continue

    def pairs_ext(self, maxhd, pairfunc = lambda node: node.isset):
        """Generator function to iterate all pairs of prefixes that are within
        a certain Hamming distance of each other.

        Yields tuples with the following format:
        (Hamming distance, prefix1, value1, prefix2, value2)
        """
        # Traverse the tree in search of pairfunc nodes.
        pairfunc_nodes = {}
        stack = [(self._root, "")]
        while stack:
            node, s = stack.pop()

            if pairfunc(node):
                pairfunc_nodes[s] = node

            for ch, next_node in node.children.iteritems():
                stack.append((next_node, s + ch))

        explored = set()
        for s, node_s in pairfunc_nodes.iteritems():
            len_s = len(s)
            explored.add((len_s, s))
            stack = [(self._root, "", 0, 0)]
            while stack:
                node, alt_s, pos, hd = stack.pop()

                if pos == len_s:
                    if pairfunc(node):
                        yield hd, s, node_s.value, \
                                alt_s, pairfunc_nodes[alt_s].value
                    # No need to go deeper into the tree because only equal
                    # length prefixes are considered.
                    continue

                n_explored = 0
                for ch, next_node in node.children.iteritems():
                    next_alt_s = alt_s + ch
                    if (len_s, next_alt_s) in explored:
                        # all pairs within maxhd for len(s) prefixes have been
                        # found already for the next_alt_s branch in the tree.
                        n_explored += 1
                        continue
                    if ch == s[pos]:
                        stack.append((next_node, next_alt_s, pos + 1, hd))
                    elif hd < maxhd:
                        stack.append((next_node, next_alt_s, pos + 1, hd + 1))

                if n_explored == len(node.children):
                    # Closing off current branch of the tree for current
                    # sequence length.
                    explored.add((len_s, alt_s))
                    continue
