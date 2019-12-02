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

        self.segmentLength = segment_length
        self.ItemType = item_dtype

        # Make a pair of empty list of segments
        self.segListPair = [ [], [] ]

        # Segment number and offset for the first available opening
        self.firstNotFullSegment = [0, 0]
        self.firstAvailableOffset = [0, 0]

        # currentSegment is used to count segments in loops over
        # segments, e.g., to push particles.
        self.currentSegment = [0, 0]

        self.nSeg = [0, 0]
        self.nPmax = [0, 0]
        for iSA in (0, 1):
            # Add the numpy array for the first segment
            self.segListPair[iSA].append(np_m.empty(self.segmentLength, dtype=item_dtype))
#            self.segListPair[iSA].append(np_m.zeros(self.segmentLength, dtype=item_dtype))
            # Count the number of segments:
            self.nSeg[iSA] = len(self.segListPair[iSA])
            # Maximum number of particles that can be stored at present
            self.nPmax[iSA] = self.nSeg[iSA]*self.segmentLength

        # Identify which of the array pairs is being written to. When
        # the particles are being initialized, we write to the first
        # one of the array pair, i.e., index 0.
        self.inSegmentedArray = 1
        self.outSegmentedArray = 0

        # Make a list of segments with erased items.  The segment items are
#        self.SegmentsWithHoles = []
#    def __init__(self, segment_length, item_dict): ENDDEF

#class SegmentedArrayPair_C(object):
    def reset_array(self):
        """Set the first available slot to the start of the segmented
           array.
           
        """
        
        self.firstAvailableSegment = 0
        self.firstAvailableOffset = 0

        return
#   def reset_array(self):ENDDEF

#class SegmentedArrayPair_C(object):
    def push_back(self, item_input):
        """Adds an item to the end of the "out" SegmentedArray.

           The item is either a tuple containing a complete structure, or a 1-element
           Numpy array. If the last segment is full, create a new segment.  This
           assumes that all the segments, except maybe the last one, are full.

           :param item_input: A tuple containing a complete item structure, or a
                              1-element Numpy array containing an item.
           :type input_item: tuple(float,...) or ndarray with 1 element.
           :var int full_index: The full index into the SA.
           :return: (offset into the "out" segment, full SA index of item)
           old: return (item structure, full SA index in the "out" array)

        """

        # Abbreviations
        outSA = self.outSegmentedArray

        # If we've reached the end of the current segment, we need to
        # switch to the next segment, if there is one, or else add a
        # new segment.
        if self.firstAvailableOffset[outSA] == self.segmentLength:
            # If another segment is already available, use
            # it. Otherwise, allocate a new segment.
            self.currentSegment[outSA] += 1
            if self.currentSegment[outSA] < self.nSeg[outSA]:
                self.firstNotFullSegment[outSA] += 1
                self.firstAvailableOffset[outSA] = 0
            else:
                # The following call increments
                # firstNotFullSegment[] and nSeg[], and sets
                # firstAvailableOffset[] = 0
                self.add_segment(outSA)

# Untested:
        if type(item_input) is np_m.ndarray:
            item = item_input[0]
        else:
            item = item_input
    
        vec = self.segListPair[outSA][self.firstNotFullSegment[outSA]]
        vec[self.firstAvailableOffset[outSA]] = item

        # Compute the full zero-based index of the particle for return
        full_index = self.firstNotFullSegment[outSA]*self.segmentLength + self.firstAvailableOffset[outSA]

        # Increment the next available slot for next time
        self.firstAvailableOffset[outSA] += 1

        return self.firstAvailableOffset[outSA]-1, full_index
# old:       return vec[self.firstAvailableOffset[outSA]-1], full_index
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

        (seg, offset) = divmod(i, self.segmentLength)
        outSA = self.outSegmentedArray
        return (self.segListPair[outSA][seg], offset)
#    def get_segment_and_offset(self, i):ENDDEF

#class SegmentedArrayPair_C(object):
    def get_item(self, i):
        """Returns a REFERENCE to the i'th item from the "out" SegmentedArray, since that's
           the up-to-date array.  The index is zero-based.

           :param int i: The full index of an item to be retrieved.

        """

        # Abbreviations
        outSA = self.outSegmentedArray

#        print 'get_item: i =', i, 'outSA = ', self.outSegmentedArray

        (seg, offset) = divmod(i, self.segmentLength)
        return self.segListPair[outSA][seg][offset]
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
        outSA = self.outSegmentedArray

        (seg, offset) = divmod(i, self.segmentLength)
        return self.segListPair[outSA][seg][offset]

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
        outSA = self.outSegmentedArray

        (seg, offset) = divmod(i, self.segmentLength)
        self.segListPair[outSA][seg][offset] = value

        return

#class SegmentedArrayPair_C(object):
    def add_segment(self, theSA):
        """Adds another Segment to the selected Segmented Array to store more items.
        """
        # Abbreviations

        self.nSeg[theSA] += 1
        self.segListPair[theSA].append(np_m.empty(self.segmentLength, dtype=self.ItemType))
#        self.segListPair[theSA].append(np_m.zeros(self.segmentLength, dtype=self.ItemType))
        self.firstNotFullSegment[theSA] += 1
        self.firstAvailableOffset[theSA] = 0

        return
    
#class SegmentedArrayPair_C(object):
    def init_out_loop(self):
        """Initialize a loop over the segments of the "out" array.

           This is used to start a loop that doesn't change the
           positions of the items, so using the "in" array isn't
           needed.

           :return: (number of items in the first segment of "out" SA,
                     ref to first segment of "out" SA)
        """

        # Abbreviations
        outSA = self.outSegmentedArray

        # This is used to count through the segments
        self.currentSegment[outSA] = 0

        # Return the first segment of the "out" array

#        print 'init_out_loop: self.firstNotFullSegment[outSA]=', self.firstNotFullSegment[outSA]

        # If the first segment is also the last segment, the number of
        # items in it may be less than the length of the segment
        if self.firstNotFullSegment[outSA] == 0:
            lastItem = self.firstAvailableOffset[outSA]
#            print 'init_out_loop: lastItem=', lastItem
            # If this is the last segment and it's empty, return None.
            # ===> This should cause the caller to break out of its loop over segments <===
            # If it's not empty, we return the non-empty items further
            # down.
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.segmentLength

        segIndex = 0
        return (lastItem, self.segListPair[outSA][segIndex][0:lastItem])
#    def init_out_loop(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def init_inout_loop(self):
        """Initialize a loop over the segments.  The loop should call
           get_next_segment("in") and get_next_out_segment() in the
           iteration loop.

           :return: (number of items in the first segment of "in" SA,
                     ref to first segment of "in" SA,
                     ref to first segment of "out" SA)
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

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
        self.currentSegment[0] = 0
        self.currentSegment[1] = 0

        # Return the first segment of the "in" array

#        print 'init_segment_loop: self.firstNotFullSegment[inSA]=', self.firstNotFullSegment[inSA]

        # If the first segment is also the last segment, the number of
        # items in it may be less than the length of the segment
        if self.firstNotFullSegment[inSA] == 0:
            lastItem = self.firstAvailableOffset[inSA]
#            print 'init_segment_loop: lastItem=', lastItem
            # If this is the last segment and it's empty, return None.
            # ===> This should cause the caller to break out of its loop over segments <===
            # If it's not empty, we return the non-empty items further
            # down.
            if lastItem == 0: return (0, None, None)
        else:
            lastItem = self.segmentLength

        # Store the member values
        self.inSegmentedArray = inSA
        self.outSegmentedArray = outSA

        segIndex = 0

        return (lastItem, self.segListPair[inSA][segIndex][0:lastItem], self.segListPair[outSA][segIndex])
#    def init_inout_loop(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_next_segment(self, InOut):
        """Returns a reference to the next segment of either the "in"
           or "out" array, and the number of active items.
        """

        if InOut == "in":
            theSA = self.inSegmentedArray
        elif InOut == "out":
            theSA = self.outSegmentedArray

        self.currentSegment[theSA] += 1
        segIndex = self.currentSegment[theSA]

        # If the segment index exceeds occupied limit, return None.
        if segIndex > self.firstNotFullSegment[theSA]: return (0, None)
        # If this is the last segment and it's empty, return None.
        # ===> This should cause the caller to break out of the loop <===
        # If it's not empty, return the non-empty items
        if self.firstNotFullSegment[theSA] == segIndex:
            lastItem = self.firstAvailableOffset[theSA]
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.segmentLength

#        self.currentSegment[inSA] += 1 # Increment for next time
        
# Q: What's this returning? Looks like a COPY of the items?  Or is it a ref?
# A: It's a ref to a numpy array with the indicated length.
#        print 'seg.py', self.segListPair[inSA][segIndex]
        return (lastItem, self.segListPair[theSA][segIndex][0:lastItem])
#    def get_next_in_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_next_out_segment(self):
        """Returns a reference to the next segment of the "out" array.
           This method looks similar to the push_back() method above, since the
           "out" array is effectively scratch space.

           If this function is called, it assumes that you need space
           to write on, so it will allocate a new segment if we're out
           of "out" segments.

           :return: reference to next segment of the "out" array.

        """

        # Abbreviations
        outSA = self.outSegmentedArray
        self.currentSegment[outSA] += 1
        segIndex = self.currentSegment[outSA]

        # If another segment is already available, use
        # it. Otherwise, allocate a new segment.
        if segIndex < self.nSeg[outSA]:
            self.firstNotFullSegment[outSA] += 1
            self.firstAvailableOffset[outSA] = 0
        else:
            # The following call increments the variables
            # firstNotFullSegment[] and nSeg[], and sets
            # firstAvailableOffset[] = 0
            self.add_segment(outSA)

        return self.segListPair[outSA][segIndex]
#    def get_next_out_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_current_out_segment(self):
        """Returns a reference to the current segment of the "out" array.

           :return: reference to current segment of the "out" array.

        """

        # Abbreviations
        outSA = self.outSegmentedArray
        segIndex = self.currentSegment[outSA]

        return self.segListPair[outSA][segIndex]
#    def get_current_out_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_number_of_segments(self):
        """Returns the current number of segments in the "in" and
           "out" arrays.
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        return (self.nSeg[inSA], self.nSeg[outSA])

#class SegmentedArrayPair_C(object):
    def get_number_of_items(self):
        """Returns the total number of items currently stored in the
           "out" array.

           :return: Number of items in the "out" SA.
        """

        # Abbreviations
        outSA = self.outSegmentedArray

        return self.firstNotFullSegment[outSA]*self.segmentLength + self.firstAvailableOffset[outSA]

#class SegmentedArrayPair_C(object):
    def set_number_of_items(self, InOut, n_items):
        """Sets the number of active items currently stored.

           :param str InOut: Either "in" or "out" depending on whether
                             we're dealing with the in or out SA.

           :return: Nothing is returned.
        """

        if InOut == "in":
            theSA = self.inSegmentedArray
        elif InOut == "out":
            theSA = self.outSegmentedArray

        # Compute the segment and offset of the last item. The
        # zero-based index of the last item is n_items-1
        (seg, offset) = divmod(n_items-1, self.segmentLength)

        self.firstNotFullSegment[theSA] = seg
        self.firstAvailableOffset[theSA] = offset+1

        return

#class SegmentedArrayPair_C(object):
    def get_capacity(self):
        """Returns the total number of items that can currently be stored.
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        return (self.nSeg[inSA]*self.segmentLength, self.nSeg[outSA]*self.segmentLength)

#class SegmentedArrayPair_C(object):
    def get_number_of_mbytes(self):
        """Returns the number of megabytes allocated for the item arrays.

           Uses the Numpy attribute nbytes.

           :return: A 2-tuple with the number of megabytes in the "in" and "out"
                    SegmentedArrays.

        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        return (self.nSeg[inSA]*self.segListPair[inSA][0].nbytes/(1.0e6),
                self.nSeg[outSA]*self.segListPair[outSA][0].nbytes/(1.0e6))

#class SegmentedArrayPair_C(object):
    def compress(self):
        """Fill in empty positions in each segment due to items that
        are no longer needed.  Should completely empty segments be
        released?
        """

        # Loop on segments
        for i in range(self.nSeg[cSA]):
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
            theSA = self.inSegmentedArray
        elif InOut == "out":
            theSA = self.outSegmentedArray

        full_index = self.currentSegment[theSA]*self.segmentLength + indx

        return full_index
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
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        full_index_in = self.currentSegment[inSA]*self.segmentLength + i_in
        full_index_out = self.currentSegment[outSA]*self.segmentLength + i_out

        return (full_index_in, full_index_out)
#    def get_full_indices(self, i_in, i_out):ENDDEF

#class SegmentedArrayPair_C(object):
    def delete_item(self, full_index):
        """Deletes an item given its full index.
        """
        # Compute the segment number
        iseg = full_index/self.segmentLength
        # Compute the offset into this segment
        ioff = full_index%self.segmentLength

        # Save the location

        self.HolesIndices[iseg][Nholes] = ioff

        vec = self.segList[self.firstNotFullSegment]

        i = self.firstAvailableIndex
        if i == segment_length:
            self.add_segment()

        i += 1
        
        pass
#    def delete_item(self, full_index): ENDDEF

#class SegmentedArrayPair_C(object): ENDCLASS
