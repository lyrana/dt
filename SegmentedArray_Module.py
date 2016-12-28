# SegmentedArray

__version__ = 0.1
__author__ = 'T. P. Hughes'
__all__ = ['SegmentedArray_C.SegList',
           'SegmentedArray_C.nSeg', ]

# use initial underscores for locals

import numpy as np_M

class SegmentedArray_C(object):
    """This class implements a segmented array of items.  Each
    segments is a NUMPY array.  Each item in these arrays is a
    struct. The fields in the struct have names and data types. All
    segments are the same length and use the same struct.

    SegList[0] is the first NUMPY array.
    SegList[n-1] is the nth NUMPY array.

    Suppose the dtype of the NUMPY array has "names" of 'x', 'y', 'z', then:

    SegList[n-1]['x'] is an array containing all the x values in the nth segment
    SegList[n-1]['x'][0] is the x value of the first item in the nth segment
    SegList[n-1][0] is a 'data record' with the the x, y, z values of the first item in the nth segment
    SegList[n-1][0]['x'] is the x value of the first item in the nth segment
    
    Print the field names of the struct in segment n: SegList[n-1].dtype.names 
    Print all the values of field 'x' in the nth segment: SegList[n-1]['x'] 
    Print the byte stride between 'x' values: SegList[0].strides
    """

#    SegmentLength = 100 # static class attribute
# Provide a classmethod to set this value?

#class SegmentedArray_C(object):
    def __init__(self, segment_length, item_dict):
        """Set up the segmented array.
        """

        self.SegmentLength = segment_length
        self.ItemDict = item_dict

        # Make an empty list of segments
        self.SegList = []

        # Add the numpy array for the first segment
        self.SegList.append(np_M.empty(self.SegmentLength, dtype=item_dict))
        
        # Segment number and offset for the first available opening
        self.FirstAvailableSegment = 0
        self.FirstAvailableOffset = 0

        # Count the number of segments:
        self.nSeg = len(self.SegList)

        # Maximum number of particles that can be stored at present
        self.nPmax = self.nSeg*self.SegmentLength

        # Make a list of segments with erased items.  The segment items are
        self.SegmentsWithHoles = []

#class SegmentedArray_C(object):
    def put(self, item_input):
        """Adds an item to the SegmentedArray, without specifying an
        index.  Creates a new Segment, if needed.
        """

        # add another segment, if needed:
        if self.FirstAvailableOffset == self.SegmentLength:
            self.add_segment()

# i += 1; self.FirstAvailableOffset += 1;  # Are both needed?
        vec = self.SegList[self.FirstAvailableSegment]
        vec[self.FirstAvailableOffset] = item_input

        # Compute the full index of the particle for return
        full_index = self.FirstAvailableSegment*self.SegmentLength + self.FirstAvailableOffset

        self.FirstAvailableOffset += 1

        return vec[self.FirstAvailableOffset-1], full_index

#class SegmentedArray_C(object):
    def get(self, i):
        """Returns a reference to the i'th item from the SegmentedArray.
        """
        (seg, offset) = divmod(i, self.SegmentLength)
        return self.SegList[seg][offset]

#class SegmentedArray_C(object):
    def __getitem__(self, i):
        """Defines subscripting, i.e., b=a[i]. Returns a reference to
        the i'th item from the SegmentedArray. Example:
          getparticle = self.seg_array_obj[i]
          print getparticle['x']
          print getparticle['weight']
        """
        (seg, offset) = divmod(i, self.SegmentLength)
        return self.SegList[seg][offset]

#class SegmentedArray_C(object):
    def __setitem__(self, i, value):
        """Defines assignments, i.e., a[i] = value. Example:
            x=0.0; y=1.0; z=2.0; ux=3.0; uy=4; uz=5.0; weight = 101.1
            putparticle = (x, y, z, ux, uy, uz, weight,)
            self.seg_array_obj[i] = putparticle
        """
        (seg, offset) = divmod(i, self.SegmentLength)
        self.SegList[seg][offset] = value
   
#class SegmentedArray_C(object):
    def add_segment(self):
        """Adds another Segment to store more items.
        """
        self.nSeg += 1
        self.SegList.append(np_M.empty(self.SegmentLength, dtype=self.ItemDict))

        # A list of the locations of unneeded items in the segment
##        self.HoleIndices[self.nSeg] = self.HoleIndices.append(np.empty(self.SegmentLength, dtype=int))
        # Initialize this list to -1 to indicate no holes.
##        self.HoleIndices[self.nSeg][:] = -1

        self.FirstAvailableSegment += 1
        self.FirstAvailableOffset = 0

#class SegmentedArray_C(object):
    def init_segment_loop(self):
        """Initialize the segment counter to 0, e.g., to start a loop that uses get_next_segment.
        """

        self.SegmentCounter = 0

#class SegmentedArray_C(object):
    def get_next_segment(self):
        """Returns a reference to the next segment of the array.
        """

        segIndex = self.SegmentCounter
        lastItem = self.SegmentLength

        # If the segment index exceeds occupied limit, return None.
        if segIndex > self.FirstAvailableSegment: return None
        # If this is the last segment and it's empty, return None.
        # ===> This should cause the caller to break out of the loop <===
        # If it's not empty, return the non-empty items
        if self.FirstAvailableSegment == segIndex:
            lastItem = self.FirstAvailableOffset
            if lastItem == 0: return None

        self.SegmentCounter += 1 # Increment for next time
        
# Q: What's this returning? Looks like a COPY of the items?  Or is it a ref?
# A: It's a ref to a numpy array with the indicated length.
#        print 'seg.py', self.SegList[segIndex]
        return self.SegList[segIndex][0:lastItem]
#    def get_next_segment(self): ENDDEF

#class SegmentedArray_C(object):
    def number_of_segments(self):
        """Returns the current number of segments.
        """
        return self.nSeg

#class SegmentedArray_C(object):
    def number_of_stored_items(self):
        """Returns the total number of items currently stored.
        """
        return self.FirstAvailableSegment*self.SegmentLength + self.FirstAvailableOffset

#class SegmentedArray_C(object):
    def item_capacity(self):
        """Returns the total number of items currently stored.
        """
        return self.nSeg*self.SegmentLength

#class SegmentedArray_C(object):
    def number_of_mbytes(self):
        """Returns the number of megabytes allocated for the item arrays.
        """
        return self.nSeg*self.SegList[0].nbytes/(1.0e6)

#class SegmentedArray_C(object):
    def compress(self):
        """Fill in empty positions in each segment due to items that
        are no longer needed.  Should completely empty segments be
        released?
        """

        # Loop on segments
        for i in range(self.nSeg):
            compress_segment(self, i)
        pass
#    def compress(self): ENDDEF

#class SegmentedArray_C(object):
    def compress_segment(self, i):
        """Fill in empty positions in a segment
        """
#    def compress_segment(self, i): ENDDEF

        
#class SegmentedArray_C(object):
    def delete_item(self, full_index):
        """Deletes an item given its full index.
        """
        # Compute the segment number
        iseg = full_index/self.SegmentLength
        # Compute the offset into this segment
        ioff = full_index%self.SegmentLength

        # Save the location

        self.HolesIndices[iseg][Nholes] = ioff

        vec = self.SegList[self.FirstAvailableSegment]

        i = self.FirstAvailableIndex
        if i == segment_length:
            self.add_segment()

        i += 1
        
        pass
#    def delete_item(self, full_index): ENDDEF

#class SegmentedArray_C(object): ENDCLASS
