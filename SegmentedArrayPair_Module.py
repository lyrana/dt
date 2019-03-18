# SegmentedArrayPair

__version__ = 0.1
__author__ = 'Copyright (C) 2016 L. D. Hughes'
__all__ = ['SegmentedArrayPair_C.SegList',
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

    SegList[0] is the first NUMPY array.
    SegList[n-1] is the nth NUMPY array.

    Suppose the dtype of the NUMPY array has "names" of 'x', 'y', 'z', then:

    SegList[n-1]['x'] is an array containing all the x values in the nth segment
    SegList[n-1]['x'][0] is the x value of the first item in the nth segment
    SegList[n-1][0] is a 'data structure' with the x, y, z values of the first item in
                    the nth segment
    SegList[n-1][0]['x'] is the x value of the first item in the nth segment
    
    Print the field names of the structure in segment n: SegList[n-1].dtype.names
    Print all the values of field 'x' in the nth segment: SegList[n-1]['x']
    Print the byte stride between 'x' values: SegList[0].strides

    """

#    SegmentLength = 100 # static class attribute
# Provide a classmethod to set this value?

#class SegmentedArrayPair_C(object):
    def __init__(self, segment_length, item_dtype):
        """Set up the segmented array.
        """

        self.SEGMENTLENGTH = segment_length
        self.ItemType = item_dtype

        # Make a pair of empty list of segments
        self.SegListPair = [ [], [] ]

        # Segment number and offset for the first available opening
        self.FirstNotFullSegment = [0, 0]
        self.FirstAvailableOffset = [0, 0]

        # CurrentSegment is used to count segments in loops over
        # segments, e.g., to push particles.
        self.CurrentSegment = [0, 0]

        self.nSeg = [0, 0]
        self.nPmax = [0, 0]
        for iSA in (0, 1):
            # Add the numpy array for the first segment
            self.SegListPair[iSA].append(np_m.empty(self.SEGMENTLENGTH, dtype=item_dtype))
#            self.SegListPair[iSA].append(np_m.zeros(self.SEGMENTLENGTH, dtype=item_dtype))
            # Count the number of segments:
            self.nSeg[iSA] = len(self.SegListPair[iSA])
            # Maximum number of particles that can be stored at present
            self.nPmax[iSA] = self.nSeg[iSA]*self.SEGMENTLENGTH

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
        
        self.FirstAvailableSegment = 0
        self.FirstAvailableOffset = 0

        return
#   def reset_array(self):ENDDEF

#class SegmentedArrayPair_C(object):
    def put(self, item_input):
        """Adds an item to the 'out' SegmentedArray, without specifying an
           index.  The item is a tuple containing a complete structure.
           Creates a new Segment, if needed.  This assumes that all
           the segments, except maybe the last one, are full.

           :param item_input: A tuple containing a complete item structure
           :type input_item: tuple(float,...)
           :var int full_index: The full index into the SA.
           :return: (item structure, full SA index)
        """

        # Abbreviations
        outSA = self.outSegmentedArray

        # If we've reached the end of the current segment, we need to
        # switch to the next segment, if there is one, or else add a
        # new segment.
        if self.FirstAvailableOffset[outSA] == self.SEGMENTLENGTH:
            # If another segment is already available, use
            # it. Otherwise, allocate a new segment.
            self.CurrentSegment[outSA] += 1
            if self.CurrentSegment[outSA] < self.nSeg[outSA]:
                self.FirstNotFullSegment[outSA] += 1
                self.FirstAvailableOffset[outSA] = 0
            else:
                # The following call increments
                # FirstNotFullSegment[] and nSeg[], and sets
                # FirstAvailableOffset[] = 0
                self.add_segment(outSA)


# i += 1; self.FirstAvailableOffset += 1;  # Are both needed?
#        print 'outSA=', outSA, 'self.FirstNotFullSegment[outSA]=', self.FirstNotFullSegment[outSA]

        vec = self.SegListPair[outSA][self.FirstNotFullSegment[outSA]]
        vec[self.FirstAvailableOffset[outSA]] = item_input

        # Compute the full zero-based index of the particle for return
        full_index = self.FirstNotFullSegment[outSA]*self.SEGMENTLENGTH + self.FirstAvailableOffset[outSA]

        # Increment the next available slot for next time
        self.FirstAvailableOffset[outSA] += 1

        return vec[self.FirstAvailableOffset[outSA]-1], full_index
#    def put(self, item_input):ENDDEF

#class SegmentedArrayPair_C(object):
    def get(self, i):
        """Returns a REFERENCE to the i'th item from the 'out'
           SegmentedArray, since that's the up-to-date array.  The
           index is zero-based.

           :param int i: The full index of an item to be retrieved.
        """

        # Abbreviations
        outSA = self.outSegmentedArray

#        print 'get: i =', i, 'outSA = ', self.outSegmentedArray

        (seg, offset) = divmod(i, self.SEGMENTLENGTH)
        return self.SegListPair[outSA][seg][offset]
#    def get(self, i):ENDDEF

#class SegmentedArrayPair_C(object):
    def __getitem__(self, i):
        """Defines subscripting on a SegmentedArrayPair_C object, i.e.,
           b=a[i]. The indexing is zero-based.

           Returns a reference to the i'th item stored in the
           segmented array.  Note that internally, it gets the item
           from the SegListPair attribute, since that's the most
           useful way to use indexing.

           Example:
             getparticle = self.seg_array_obj[i]
           So now you can write:
             print getparticle['x']
             print getparticle['weight']
        """
        # Abbreviations
        # By default, __getitem__() uses the 'out' array since it contains the most up-to-date values
        outSA = self.outSegmentedArray

        (seg, offset) = divmod(i, self.SEGMENTLENGTH)
        return self.SegListPair[outSA][seg][offset]

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
        # By definition, __setitem__() uses the 'out' array
        outSA = self.outSegmentedArray

        (seg, offset) = divmod(i, self.SEGMENTLENGTH)
        self.SegListPair[outSA][seg][offset] = value

        return

#class SegmentedArrayPair_C(object):
    def add_segment(self, theSA):
        """Adds another Segment to the selected Segmented Array to store more items.
        """
        # Abbreviations

        self.nSeg[theSA] += 1
        self.SegListPair[theSA].append(np_m.empty(self.SEGMENTLENGTH, dtype=self.ItemType))
#        self.SegListPair[theSA].append(np_m.zeros(self.SEGMENTLENGTH, dtype=self.ItemType))
        self.FirstNotFullSegment[theSA] += 1
        self.FirstAvailableOffset[theSA] = 0

        return
    
#class SegmentedArrayPair_C(object):
    def init_out_loop(self):
        """Initialize a loop over the segments of the 'out' array.

           This is used to start a loop that doesn't change the
           positions of the items, so using the 'in' array isn't
           needed.
        """

        # Abbreviations
        outSA = self.outSegmentedArray

        # These are used to count through the segments
        self.CurrentSegment[outSA] = 0

        # Return the first segment of the 'out' array

#        print 'init_out_loop: self.FirstNotFullSegment[outSA]=', self.FirstNotFullSegment[outSA]

        # If the first segment is also the last segment, the number of
        # items in it may be less than the length of the segment
        if self.FirstNotFullSegment[outSA] == 0:
            lastItem = self.FirstAvailableOffset[outSA]
#            print 'init_out_loop: lastItem=', lastItem
            # If this is the last segment and it's empty, return None.
            # ===> This should cause the caller to break out of its loop over segments <===
            # If it's not empty, we return the non-empty items further
            # down.
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.SEGMENTLENGTH

        segIndex = 0
        return (lastItem, self.SegListPair[outSA][segIndex][0:lastItem])
#    def init_out_loop(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def init_inout_loop(self):
        """Initialize a loop over the segments.  The loop should use
           get_next_segment('in') and get_next_out_segment() in the
           iteration loop.

           :return: (number of items in the first segment of 'in' SA,
                     ref to first segment of 'in' SA)
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

#        print 'init_segment_loop A: inSA =', inSA

        # Swap the two particle storage arrays, so that the current
        # 'out' array becomes the 'in' array.
        # The 'in' array appears on the RHS of expressions, i.e., is
        # input to an expression.
        # The 'out' array appears on the LHS of expressions, i.e., is
        # the result of an expression. It's like a scratch array: it's
        # initial values don't matter, just it's length.

        inSA = (inSA+1) % 2
        outSA = (outSA+1) % 2

#        print 'init_segment_loop B: inSA =', inSA, 'outSA =', outSA

        # These counters are used to count through the segments.
        # Segment indexing is zero-based.
        self.CurrentSegment[0] = 0
        self.CurrentSegment[1] = 0

        # Return the first segment of the 'in' array

#        print 'init_segment_loop: self.FirstNotFullSegment[inSA]=', self.FirstNotFullSegment[inSA]

        # If the first segment is also the last segment, the number of
        # items in it may be less than the length of the segment
        if self.FirstNotFullSegment[inSA] == 0:
            lastItem = self.FirstAvailableOffset[inSA]
#            print 'init_segment_loop: lastItem=', lastItem
            # If this is the last segment and it's empty, return None.
            # ===> This should cause the caller to break out of its loop over segments <===
            # If it's not empty, we return the non-empty items further
            # down.
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.SEGMENTLENGTH

        # Store the member values
        self.inSegmentedArray = inSA
        self.outSegmentedArray = outSA

        segIndex = 0

        return (lastItem, self.SegListPair[inSA][segIndex][0:lastItem], self.SegListPair[outSA][segIndex])

#    def init_inout_loop(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_next_segment(self, InOut):
        """Returns a reference to the next segment of either the 'in'
           or 'out' array, and the number of active items.
        """

        if InOut == 'in':
            theSA = self.inSegmentedArray
        elif InOut == 'out':
            theSA = self.outSegmentedArray

        self.CurrentSegment[theSA] += 1
        segIndex = self.CurrentSegment[theSA]

        # If the segment index exceeds occupied limit, return None.
        if segIndex > self.FirstNotFullSegment[theSA]: return (0, None)
        # If this is the last segment and it's empty, return None.
        # ===> This should cause the caller to break out of the loop <===
        # If it's not empty, return the non-empty items
        if self.FirstNotFullSegment[theSA] == segIndex:
            lastItem = self.FirstAvailableOffset[theSA]
            if lastItem == 0: return (0, None)
        else:
            lastItem = self.SEGMENTLENGTH

#        self.CurrentSegment[inSA] += 1 # Increment for next time
        
# Q: What's this returning? Looks like a COPY of the items?  Or is it a ref?
# A: It's a ref to a numpy array with the indicated length.
#        print 'seg.py', self.SegListPair[inSA][segIndex]
        return (lastItem, self.SegListPair[theSA][segIndex][0:lastItem])
#    def get_next_in_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_next_out_segment(self):
        """Returns a reference to the next segment of the 'out' array.
           This method is similar to the put() method above, since the
           'out' array is effectively scratch space.

           :return: reference to next segment of the 'out' array.
        """

        # Abbreviations
        outSA = self.outSegmentedArray
        self.CurrentSegment[outSA] += 1
        segIndex = self.CurrentSegment[outSA]

        # If another segment is already available, use
        # it. Otherwise, allocate a new segment.
        if segIndex < self.nSeg[outSA]:
            self.FirstNotFullSegment[outSA] += 1
            self.FirstAvailableOffset[outSA] = 0
        else:
            # The following call increments the variables
            # FirstNotFullSegment[] and nSeg[], and sets
            # FirstAvailableOffset[] = 0
            self.add_segment(outSA)

        return self.SegListPair[outSA][segIndex]
#    def get_next_out_segment(self): ENDDEF

#class SegmentedArrayPair_C(object):
    def get_number_of_segments(self):
        """Returns the current number of segments in the 'in' and
           'out' arrays.
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        return (self.nSeg[inSA], self.nSeg[outSA])

#class SegmentedArrayPair_C(object):
    def get_number_of_items(self):
        """Returns the total number of items currently stored in the
           'out' array.

           :return: Number of items in the 'out' SA.
        """

        # Abbreviations
        outSA = self.outSegmentedArray

        return self.FirstNotFullSegment[outSA]*self.SEGMENTLENGTH + self.FirstAvailableOffset[outSA]

#class SegmentedArrayPair_C(object):
    def set_number_of_items(self, InOut, n_items):
        """Sets the number of active items currently stored.

           :param str InOut: Either 'in' or 'out' depending on whether
                             we're dealing with the in or out SA.

           :return: Nothing is returned.
        """

        if InOut == 'in':
            theSA = self.inSegmentedArray
        elif InOut == 'out':
            theSA = self.outSegmentedArray

        # Compute the segment and offset of the last item. The
        # zero-based index of the last item is n_items-1
        (seg, offset) = divmod(n_items-1, self.SEGMENTLENGTH)

        self.FirstNotFullSegment[theSA] = seg
        self.FirstAvailableOffset[theSA] = offset+1

        return

#class SegmentedArrayPair_C(object):
    def get_item_capacity(self):
        """Returns the total number of items currently stored.
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        return (self.nSeg[inSA]*self.SEGMENTLENGTH, self.nSeg[outSA]*self.SEGMENTLENGTH)
#        return self.nSeg[cSA]*self.SEGMENTLENGTH

#class SegmentedArrayPair_C(object):
    def get_number_of_mbytes(self):
        """Returns the number of megabytes allocated for the item arrays.
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        return (self.nSeg[inSA]*self.SegListPair[inSA][0].nbytes/(1.0e6),
                self.nSeg[outSA]*self.SegListPair[outSA][0].nbytes/(1.0e6))

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
    def get_full_indices(self, i_in, i_out):
        """
           This computes the full indices of a particular item in the
           'in' and 'out' SAs.

           :param int i_in: offset of an item into the 'in' SA
           :param int i_out: offset of an item into the 'out' SA

           :return: (full index in 'in' array, full index in 'out' array)
        """

        # Abbreviations
        inSA = self.inSegmentedArray
        outSA = self.outSegmentedArray

        full_index_in = self.CurrentSegment[inSA]*self.SEGMENTLENGTH + i_in
        full_index_out = self.CurrentSegment[outSA]*self.SEGMENTLENGTH + i_out

        return (full_index_in, full_index_out)
#    def get_full_indices(self, i_in, i_out):ENDDEF

#class SegmentedArrayPair_C(object):
    def get_full_index(self, indx, InOut):
        """
           This computes the full index of a particular item in either the 'in' or 'out'
           SAs.

           :param int indx: offset of an item into the SA
           :param str InOut: 'in' or 'out', depending of which SA is intended

           :return: full index in 'in' or 'out' array
        """

        if InOut == 'in':
            theSA = self.inSegmentedArray
        elif InOut == 'out':
            theSA = self.outSegmentedArray

        full_index = self.CurrentSegment[theSA]*self.SEGMENTLENGTH + indx

        return full_index
#    def get_full_index(self, indx, inout):ENDDEF
        
#class SegmentedArrayPair_C(object):
    def delete_item(self, full_index):
        """Deletes an item given its full index.
        """
        # Compute the segment number
        iseg = full_index/self.SEGMENTLENGTH
        # Compute the offset into this segment
        ioff = full_index%self.SEGMENTLENGTH

        # Save the location

        self.HolesIndices[iseg][Nholes] = ioff

        vec = self.SegList[self.FirstNotFullSegment]

        i = self.FirstAvailableIndex
        if i == segment_length:
            self.add_segment()

        i += 1
        
        pass
#    def delete_item(self, full_index): ENDDEF

#class SegmentedArrayPair_C(object): ENDCLASS
