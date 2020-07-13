# SegmentedArrayPair

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['SegmentedArrayPair_C.segList',
           'SegmentedArrayPair_C.nSeg', ]

# use initial underscores for locals

import sys
import numpy as np_m

class SegmentedArrayPair_C(object):
    """This class implements a segmented array of items.  Each
    segment is a Numpy structured array, i.e., each item in the array is a
    structure. The data fields in the structure have names and data types. All
    segments are the same length and use the same structure type.

    The segments are always full of structures, except possibly the
    last one.

    segList[0] is the first NUMPY array.
    segList[n-1] is the nth NUMPY array.

    Suppose the dtype of the NUMPY array has "names" of 'x', 'y', 'z', then:

    segList[n-1]['x'] is an array containing all the x values in the nth segment
    segList[n-1]['x'][0] is the x value of the first item in the nth segment
    segList[n-1][0] is a 'data structure' with the x, y, z values of the first item in
                    the nth segment
    segList[n-1][0]['x'] is the x value of the first item in the nth segment
    
    Print the field names of the structure in segment n: segList[n-1].dtype.names
    Print all the values of field 'x' in the nth segment: segList[n-1]['x']
    Print the byte stride between 'x' values: segList[0].strides

    """

#    SegmentLength = 100 # static class attribute
# Provide a classmethod to set this value?

#class SegmentedArrayPair_C(object):
    def __init__(self, segment_length, item_dtype):
        """Set up the segmented array.
        """

        self.segment_length = segment_length
        self.item_dtype = item_dtype

        # Make a pair of empty lists of segments
        self.seg_list_pair = [ [], [] ]

        # Segment number and offset for the first available opening
        #     self.first_not_full_segment() marks the last existing segment that has
        #     active particles. It's used by push_back() to add more particles.
        self.first_not_full_segment = [0, 0]
        self.first_available_offset = [0, 0]

        # currentSegment is used to count segments in loops over
        # segments, e.g., to push particles.
        self.current_segment = [0, 0]

        self.nseg = [0, 0]
        self.npmax = [0, 0]
        for iSA in (0, 1):
            # Add the numpy array for the first segment
            self.seg_list_pair[iSA].append(np_m.empty(self.segment_length, dtype=item_dtype))
#            self.seg_list_pair[iSA].append(np_m.zeros(self.segment_length, dtype=item_dtype))
            # Count the number of segments:
            self.nseg[iSA] = len(self.seg_list_pair[iSA])
            # Maximum number of particles that can be stored at present
            self.npmax[iSA] = self.nseg[iSA]*self.segment_length

        # Identify which of the array pairs is being written to. When
        # the particles are being initialized, we write to the first
        # one of the array pair, i.e., index 0.
        self.in_segmented_array = 1
        self.out_segmented_array = 0

        # Make a list of segments with erased items.  The segment items are
#        self.SegmentsWithHoles = []
#    def __init__(self, segment_length, item_dict): ENDDEF

#class SegmentedArrayPair_C(object):
    def push_back(self, item_input):
        """Adds one item to the end of the "out" SegmentedArray.

           The item is either a tuple containing a complete item structure, or a 1-element
           Numpy array that contains a complete item structure.

           If the last segment is full, create a new segment.  This assumes that all the
           segments, except maybe the last one, are full.

           :param item_input: A tuple containing a complete item structure, or a
                              1-element Numpy array containing a complete item.
           :type input_item: tuple(float,...) or ndarray with 1 element.
           :var int fullIndex: The full index into the SA.
           :return: (offset into the "out" segment, full SA index of item)
           old: return (item structure, full SA index in the "out" array)

        """

        # Abbreviations
        outSA = self.out_segmented_array

        # If we've reached the end of the current segment, we need to switch to the next
        # segment, if there is one, or else add a new segment.
        if self.first_available_offset[outSA] == self.segment_length:
            # If another segment is already available, use it. Otherwise, allocate a new
            # segment.
            self.current_segment[outSA] += 1
            if self.current_segment[outSA] < self.nseg[outSA]:
                self.first_not_full_segment[outSA] += 1
                self.first_available_offset[outSA] = 0
            else:
                # The following call increments firstNotFullSegment[] and nSeg[], and sets
                # firstAvailableOffset[] = 0
                self.add_segment(outSA)
            # These two counters should be the same:                
            assert (self.first_not_full_segment[outSA] == self.current_segment[outSA]), "first_not_full_segment = %d, current_segment = %d" % (self.first_not_full_segment[outSA], self.current_segment[outSA])
                
# Untested:
        if type(item_input) is np_m.ndarray:
            item = item_input[0]
        else:
            item = item_input
    
        vec = self.seg_list_pair[outSA][self.first_not_full_segment[outSA]]
        vec[self.first_available_offset[outSA]] = item

        # Compute the full zero-based index of the particle for return
        fullIndex = self.first_not_full_segment[outSA]*self.segment_length + self.first_available_offset[outSA]

        # Increment the next available slot for next time
        self.first_available_offset[outSA] += 1

        return self.first_available_offset[outSA]-1, fullIndex
# old:       return vec[self.first_available_offset[outSA]-1], fullIndex
#    def push_back(self, item_input):ENDDEF

#class SegmentedArrayPair_C(object):
    def get_segment_and_offset(self, i):
        """Given the full "out" SegmentedArray index of a stored item, return the Numpy array
           (segment) and the offset of the item in the "out" array.  This is used,
           e.g., to provide access to the item in Python.

           Using this instead of get_item() (see below) gives the Python caller a Numpy
           structured item instead of a dict.

           :param int i: The full index of an item to be retrieved.
           :return: The tuple (Numpy array, item offset)

        """

        (seg, offset) = divmod(i, self.segment_length)
        outSA = self.out_segmented_array
        return (self.seg_list_pair[outSA][seg], offset)
#    def get_segment_and_offset(self, i):ENDDEF

#class SegmentedArrayPair_C(object):
    def get_item(self, i):
        """Returns a REFERENCE to the i'th item from the "out" SegmentedArray, since that's
           the up-to-date array.  The index is zero-based.

           :param int i: The full index of an item to be retrieved.

        """

        # Abbreviations
        outSA = self.out_segmented_array

#        print 'get_item: i =', i, 'outSA = ', self.out_segmented_array

        (seg, offset) = divmod(i, self.segment_length)
        return self.seg_list_pair[outSA][seg][offset]
#    def get_item(self, i):ENDDEF

#class SegmentedArrayPair_C(object):
    def __getitem__(self, i):
        """Defines subscripting on a SegmentedArrayPair_C object, i.e.,
           b=a[i]. The indexing is zero-based.

           Returns a reference to the i'th item stored in the
           segmented array.  Note that internally, it gets the item
           from the segListPair attribute, since that's the most
           useful way to use indexing.

           Example:
             getparticle = self.seg_array_pair_obj[i]
           So now you can write:
             print getparticle['x']
             print getparticle['weight']
        """
        # Abbreviations
        # By default, __getitem__() uses the "out" array since it contains the most up-to-date values
        outSA = self.out_segmented_array

        (seg, offset) = divmod(i, self.segment_length)
        return self.seg_list_pair[outSA][seg][offset]

#class SegmentedArrayPair_C(object):
    def __setitem__(self, i, value):
        """Defines assignments, i.e., a[i] = value.
           Assumes that storage location 'i' already exists.
           Example:
             x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
             putparticle = (x, y, z, ux, uy, uz, weight,)
             self.seg_array_obj[i] = putparticle
        """

        # Abbreviations
        # By definition, __setitem__() uses the "out" array
        outSA = self.out_segmented_array

        (seg, offset) = divmod(i, self.segment_length)
        self.seg_list_pair[outSA][seg][offset] = value

        return

#class SegmentedArrayPair_C(object):
    def add_segment(self, theSA):
        """Adds another Segment to the selected Segmented Array to store more items.
        """
        # Abbreviations
        # none

        # Add a new segment to the Segmented Array
        self.nseg[theSA] += 1
        self.seg_list_pair[theSA].append(np_m.empty(self.segment_length, dtype=self.item_dtype))
#        self.seg_list_pair[theSA].append(np_m.zeros(self.segment_length, dtype=self.item_dtype))

        # Update the capacity of the Segmented Array
        self.npmax[theSA] = self.nseg[theSA]*self.segment_length

        self.first_not_full_segment[theSA] += 1
        self.first_available_offset[theSA] = 0

        return
    
#class SegmentedArrayPair_C(object):
    def init_out_loop(self):
        """Initialize a loop over the segments of the "out" array.

           This is used to start a loop over the "out" array that doesn't change the "out"
           array.

           :return: (number of items in the first segment of "out" SA,
                     ref to first segment of "out" SA)

        """

        # Abbreviations
        outSA = self.out_segmented_array

        # This is used to count through the segments
        self.current_segment[outSA] = 0

        # Return the first segment of the "out" array

#        print 'init_out_loop: self.first_not_full_segment[outSA]=', self.first_not_full_segment[outSA]

        # If the first segment is also the last segment, the number of
        # items in it may be less than the length of the segment
        if self.first_not_full_segment[outSA] == 0:
            lastItem = self.first_available_offset[outSA]
#            print 'init_out_loop: lastItem=', lastItem
            # If this is the last segment and it's empty, return None.
            # ===> This should cause the caller to break out of its loop over segments <===
            # If it's not empty, we return the non-empty items further
            # down.
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.segment_length

        segIndex = 0
        return (lastItem, self.seg_list_pair[outSA][segIndex][0:lastItem])
#    def init_out_loop(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def init_inout_loop(self):
        """Initialize a loop that reads from the "in" segments and writes to the "out" segments.

           The loop should call get_next_segment("in") and get_next_out_segment() in the
           iteration loop.

           The loop resets self.first_not_full_segment[outSA] to zero for the "out"
           segment.

           :return: (number of items in the first segment of "in" SA,
                     ref to first segment of "in" SA,
                     ref to first segment of "out" SA)

        """

        # Abbreviations
        inSA = self.in_segmented_array
        outSA = self.out_segmented_array

#        print 'init_segment_loop A: inSA =', inSA

        # Swap the two particle storage arrays, so that the current
        # "out" array becomes the "in" array.
        # The "in" array appears on the RHS of expressions, i.e., is
        # input to an expression.
        # The "out" array appears on the LHS of expressions, i.e., is
        # the result of an expression. It's like a scratch array: it's
        # initial values don't matter, just it's length.

        inSA = (inSA+1) % 2
        outSA = (outSA+1) % 2

#        print 'init_segment_loop B: inSA =', inSA, 'outSA =', outSA

        # These counters are used to count through the segments.
        # Segment indexing is zero-based.
        self.current_segment[0] = 0
        self.current_segment[1] = 0

        # These give next available storage location in the "out" array.
        self.first_not_full_segment[outSA] = 0
        self.first_available_offset[outSA] = 0
        
        # Return the first segment of the "in" array

#        print 'init_segment_loop: self.first_not_full_segment[inSA]=', self.first_not_full_segment[inSA]

        # If the first segment is also the last segment, the number of
        # items in it may be less than the length of the segment
        if self.first_not_full_segment[inSA] == 0:
            lastItem = self.first_available_offset[inSA]
#            print 'init_segment_loop: lastItem=', lastItem
            # If this is the last segment and it's empty, return None.
            # ===> This should cause the caller to break out of its loop over segments <===
            # If it's not empty, we return the non-empty items further
            # down.
            if lastItem == 0: return (0, None, None)
        else:
            lastItem = self.segment_length

        # Store the member values
        self.in_segmented_array = inSA
        self.out_segmented_array = outSA

        segIndex = 0

        return (lastItem, self.seg_list_pair[inSA][segIndex][0:lastItem], self.seg_list_pair[outSA][segIndex])
#    def init_inout_loop(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_next_segment(self, InOut):
        """Returns the number of active items in the next segment of either the "in" or "out"
           array, and a reference to the items.

           The variable self.first_not_full_segment[theSA] is the zero-based index of the
           last segment that has active particles. This variable is not changed.

           Unlike get_next_out_segment(), this function does not allocate new segments. If
           there's no next segment, it returns a zero for the number of items.

        """

        if InOut == "in":
            theSA = self.in_segmented_array
        elif InOut == "out":
            theSA = self.out_segmented_array

        self.current_segment[theSA] += 1
        segIndex = self.current_segment[theSA]

        # If the segment index exceeds occupied limit, return None.
        if segIndex > self.first_not_full_segment[theSA]: return (0, None)
        # If this is the last segment and it's empty, return None.
        # ===> This should cause the caller to break out of the loop <===
        # If it's not empty, return the non-empty items
        if self.first_not_full_segment[theSA] == segIndex:
            lastItem = self.first_available_offset[theSA]
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.segment_length

#        self.current_segment[inSA] += 1 # Increment for next time
        
# Q: What's this returning? Looks like a COPY of the items?  Or is it a ref?
# A: It's a ref to a numpy array with the indicated length.
#        print 'seg.py', self.seg_list_pair[inSA][segIndex]
        return (lastItem, self.seg_list_pair[theSA][segIndex][0:lastItem])
#    def get_next_in_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_next_out_segment(self):
        """Returns a reference to the next segment of the "out" array.

           If this function is called, it assumes that you need space to write on, so it
           will allocate a new segment if we're out of "out" segments.

           The variable self.first_not_full_segment[outSA] is incremented to give the
           zero-based index of the returned segment. This allows push_back() to be called
           to add new particles.

           This method looks similar to the push_back() method above, since the "out"
           array is effectively scratch space that can be written to.

           :return: reference to next segment of the "out" array.

        """

        # Abbreviations
        outSA = self.out_segmented_array
        self.current_segment[outSA] += 1
        segIndex = self.current_segment[outSA]

        # If another segment is already available, use it. Otherwise, allocate a new
        # segment.
        if segIndex < self.nseg[outSA]:
            # Increment self.first_not_full_segment[outSA] if the new value is the index
            # of an existing segment.
            if self.first_not_full_segment[outSA] + 1 < self.nseg[outSA]:
                self.first_not_full_segment[outSA] += 1
                self.first_available_offset[outSA] = 0
        else:
            # The following call increments the variables firstNotFullSegment[] and
            # nSeg[], and sets firstAvailableOffset[] = 0
            self.add_segment(outSA)
        # These two counters should be the same:
        assert (self.first_not_full_segment[outSA] == segIndex), "first_not_full_segment = %d, segIndex = %d" % (self.first_not_full_segment[outSA], segIndex)

        return self.seg_list_pair[outSA][segIndex]
#    def get_next_out_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_current_out_segment(self):
        """Returns a reference to the current segment of the "out" array.

           :return: reference to current segment of the "out" array.

        """

        # Abbreviations
        outSA = self.out_segmented_array
        segIndex = self.current_segment[outSA]

        return self.seg_list_pair[outSA][segIndex]
#    def get_current_out_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_number_of_segments(self):
        """Returns the current number of segments in the "in" and
           "out" arrays.
        """

        # Abbreviations
        inSA = self.in_segmented_array
        outSA = self.out_segmented_array

        return (self.nseg[inSA], self.nseg[outSA])

#class SegmentedArrayPair_C(object):
    def get_number_of_items(self):
        """Returns the total number of items currently stored in the
           "out" array.

           :return: Number of items in the "out" SA.
        """

        # Abbreviations
        outSA = self.out_segmented_array

        return self.first_not_full_segment[outSA]*self.segment_length + self.first_available_offset[outSA]

#class SegmentedArrayPair_C(object):
    def set_number_of_items(self, InOut, n_items):
        """Sets the number of active items currently stored.

           :param str InOut: Either "in" or "out" depending on whether
                             we're dealing with the in or out SA.

           :return: Nothing is returned.
        """

        if InOut == "in":
            theSA = self.in_segmented_array
        elif InOut == "out":
            theSA = self.out_segmented_array

        # Compute the segment and offset of the last item. The
        # zero-based index of the last item is n_items-1
        (seg, offset) = divmod(n_items-1, self.segment_length)

        self.first_not_full_segment[theSA] = seg
        self.first_available_offset[theSA] = offset+1

        return

#class SegmentedArrayPair_C(object):
    def get_capacity(self):
        """Returns the total number of items that can currently be stored.
        """

        # Abbreviations
        inSA = self.in_segmented_array
        outSA = self.out_segmented_array

        return (self.nseg[inSA]*self.segment_length, self.nseg[outSA]*self.segment_length)

#class SegmentedArrayPair_C(object):
    def get_number_of_mbytes(self):
        """Returns the number of megabytes allocated for the item arrays.

           Uses the Numpy attribute nbytes.

           :return: A 2-tuple with the number of megabytes in the "in" and "out"
                    SegmentedArrays.

        """

        # Abbreviations
        inSA = self.in_segmented_array
        outSA = self.out_segmented_array

        return (self.nseg[inSA]*self.seg_list_pair[inSA][0].nbytes/(1.0e6),
                self.nseg[outSA]*self.seg_list_pair[outSA][0].nbytes/(1.0e6))

#class SegmentedArrayPair_C(object):
    def compress(self):
        """Fill in empty positions in each segment due to items that
        are no longer needed.  Should completely empty segments be
        released?
        """

        # Loop on segments
        for i in range(self.nseg[cSA]):
            compress_segment(self, i)
        pass
#    def compress(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def compress_segment(self, i):
        """Fill in empty positions in a segment
        """
#    def compress_segment(self, i): ENDDEF


#class SegmentedArrayPair_C(object):
    def get_full_index(self, indx, InOut):
        """
           This computes the full index of a particular item in either the "in" or "out"
           SAs.

           :param int indx: offset of an item into the SA
           :param str InOut: "in" or "out", depending of which SA is intended

           :return: full index in "in" or "out" array
        """

        if InOut == "in":
            theSA = self.in_segmented_array
        elif InOut == "out":
            theSA = self.out_segmented_array

        fullIndex = self.current_segment[theSA]*self.segment_length + indx

        return fullIndex
#    def get_full_index(self, indx, inout):ENDDEF
        
#class SegmentedArrayPair_C(object):
    def get_full_indices(self, i_in, i_out):
        """
           This computes the full indices of a particular item in the
           "in" and "out" SAs.

           :param int i_in: offset of an item into the "in" SA
           :param int i_out: offset of an item into the "out" SA

           :return: (full index in "in" array, full index in "out" array)
        """

        # Abbreviations
        inSA = self.in_segmented_array
        outSA = self.out_segmented_array

        fullIndexIn = self.current_segment[inSA]*self.segment_length + i_in
        fullIndexInOut = self.current_segment[outSA]*self.segment_length + i_out

        return (fullIndexIn, fullIndexInOut)
#    def get_full_indices(self, i_in, i_out):ENDDEF


#class SegmentedArrayPair_C(object): ENDCLASS
