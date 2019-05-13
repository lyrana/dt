/*! \file SegmentedArrayPair.h

  \brief This file has the source code for a C++ implementation of SegmentedArray storage.

  \namespace dnt

*/

#ifndef SEGMENTEDARRAYPAIR_H
#define SEGMENTEDARRAYPAIR_H

#include <iostream>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include "Pstruct.h"

//using namespace std; 

//template<typename PS>
//using SegList = std::vector<py::array_t<PS,0>>;

namespace dnt
{
  /*! \class WhichArray

      \brief WhichArray is an enum class to tag the 'in' and 'out' arrays of a
             SegmentedArrayPair

  */
  enum class WhichArray
  {
    in_array,
    out_array
  };
  
/*! \class SegmentedArrayPair
    \brief The SegmentedArrayPair class is a C++ version of the Python class SegmentedArrayPair_C.

    This class implements a segmented array of items.  Each segment is
    a Numpy structured array, i.e., each item in the array is a
    struct. The data fields in the struct have names and data
    types. All segments are the same length and use the same struct
    type.

    The class is templated on the struct (PT) that it has to
    store.

    The segments are always full of structs, except possibly the
    last one.

    SegList[0] is the first Numpy array.
    SegList[n-1] is the nth Numpy array.

    \param PT is the type of particle struct.
    \sa Ptype, Pstruct, WhichArray

*/
  template<Ptype PT>
  class SegmentedArrayPair
  {

  private:
//    int segmentLength;
    py::ssize_t segmentLength;

    // Make a pair of empty list of segments

    // This contains two lists of Numpy arrays
//    SegList segListPair[2];
    std::vector<py::array_t<Pstruct<PT>,0>> segListPair[2];

    unsigned nSeg[2];
    unsigned nPmax[2];
    
    // Segment number and offset for the first available opening
    unsigned firstNotFullSegment[2];
    unsigned firstAvailableOffset[2];
    
    // currentSegment is used to count segments in loops over
    // segments, e.g., to push particles.
    unsigned currentSegment[2];

    // Identify which of the array pairs is being written to. When
    // the particles are being initialized, we write to the first
    // one of the array pair, i.e., index 0.
    unsigned inSegmentedArray;
    unsigned outSegmentedArray;


    //! Create a Numpy array in C++.
    /*!

      allocate_a_segment() creates a Numpy array of length segmentLength for the PT
      particle type, and returns a pointer to it.

      \param segmentLength is the length of the Numpy array to be allocated
      \return pointer to the new Numpy array

    */
    py::array_t<Pstruct<PT>, 0> allocate_a_segment(const py::ssize_t segmentLength) 
//    py::array_t<PS> allocate_a_segment(const int segmentLength) 
    {
      // No pointer is passed, so Numpy will allocate the buffer
      return py::array_t<Pstruct<PT>, 0>(segmentLength);
    }

    //! Add another segment to the specified SegmentedArray in order to store more items.
    /*!

      \param theSA (0 or 1) identifies which of the two SegmentedArrays needs a new segment.
      \return void

    */
    void add_segment(unsigned theSA)
    {
      nSeg[theSA] += 1;
      py::array_t<Pstruct<PT>, 0> newSegment = allocate_a_segment(segmentLength);
      segListPair[theSA].push_back(newSegment);

      firstNotFullSegment[theSA] += 1;
      firstAvailableOffset[theSA] = 0;
    }
    
  public:

    //! The one and only ctor.
    /*!

      \param segment_length is the number of structs stored in one segment of the array.

      \sa Ptype

    */
    SegmentedArrayPair(py::ssize_t segment_length):
      segmentLength(segment_length),
      nSeg{0, 0},
      nPmax{0, 0},
      firstNotFullSegment{0, 0},
      firstAvailableOffset{0, 0},
      currentSegment{0, 0},
      inSegmentedArray(1),
      outSegmentedArray(0)
      {
        std::cout << "Hello from the SegmentedArrayPair ctor" << std::endl;
        for (auto iSA : {0, 1})
          {
            // Add the Numpy arrays for the first segment pair
            //            SegListPair[iSA].append(np_m.empty(self.SEGMENTLENGTH, dtype=item_dtype))
        
            //          seg_list_pair[iSA] = new SegList;
            //          seg_list_pair[iSA] = new std::vector<py::array_t<PS,0>>;
            // Create the structured array
            //          std::cout << "ctor segmentLength is " << segmentLength << std::endl;
            py::array_t<Pstruct<PT>, 0> newSegment = allocate_a_segment(segmentLength);
            segListPair[iSA].push_back(newSegment);

            //          std::cout << "type of the item in the vector is"
        
            //            self.SegListPair[iSA].append(np_m.zeros(self.SEGMENTLENGTH, dtype=item_dtype))
            // Count the number of segments:
            nSeg[iSA] = segListPair[iSA].size();
            // Maximum number of particles that can be stored at present
            nPmax[iSA] = nSeg[iSA]*segmentLength;
          }
      }

    // The dtor (See pybind 8.5 Non-public destructors)
    ~SegmentedArrayPair()
      {
        std::cout << "The dtor ~SegmentedArrayPair has been called" << std::endl;
      
        // Release the Numpy arrays?  No: this causes a crash. The
        // syntax may be wrong below. How should this memory be released?
      
        /* for (auto iSA : {0, 1}) */
        /* { */
        /*   for (auto const& seg : segListPair[iSA]) */
        /*   { */
        /*     delete &seg; // Is this correct? */
        /*   }       */
        /* } */
      } 

    //! Return the segment length
    py::ssize_t get_segment_length()
      {
        return segmentLength;
      }
    
    //! Return a py::tuple with the current number of segments in the 'in' and 'out' arrays.
    py::tuple get_number_of_segments()
      {
        // Abbreviations
        auto inSA = inSegmentedArray;
        auto outSA = outSegmentedArray;
        
        return py::make_tuple(nSeg[inSA], nSeg[outSA]);
      }


    //! Return the number of items stored in the 'out' array
      py::ssize_t get_number_of_items()
      {
        // Abbreviations
        auto outSA = outSegmentedArray;

        return firstNotFullSegment[outSA]*segmentLength + firstAvailableOffset[outSA];
      }

    //! Return the total number of items that can currently be stored.
    py::tuple get_capacity()
      {
        // Abbreviations
        auto inSA = inSegmentedArray;
        auto outSA = outSegmentedArray;

        return py::make_tuple(nSeg[inSA]*segmentLength, nSeg[outSA]*segmentLength);
      }

    //! Return the number of megabytes allocated for the item arrays
    py::tuple get_number_of_mbytes()
      {
        // Abbreviations
        auto inSA = inSegmentedArray;
        auto outSA = outSegmentedArray;

        // Get the number of bytes in one segment (a structured Numpy array)
        py::buffer_info pseg_info = segListPair[outSA][0].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        py::ssize_t itembytes = pseg_info.itemsize;
        py::ssize_t nitems = pseg_info.size;
        py::ssize_t nbytes = nitems*itembytes;
        
        return py::make_tuple(nSeg[inSA]*nbytes/(1.0e6), nSeg[outSA]*nbytes/(1.0e6));
      }
    
    //! Add an item to the 'out' SegmentedArray.
    /*!

      \param item_input is a tuple containing a complete item structure
      \return the tuple: (array member with the new item, full index)

    */
    py::tuple push_back(py::tuple item_input)
      {
        // Abbreviations
        auto outSA = outSegmentedArray;

        /*        
                  If we've reached the end of the current segment, we need to
                  switch to the next segment, if there is one, or else add a new
                  segment.
        */
        if (firstAvailableOffset[outSA] == segmentLength)
          {
            // If another segment is already available, use
            // it. Otherwise, allocate a new segment.
            currentSegment[outSA] += 1;
            if (currentSegment[outSA] < nSeg[outSA])
              {
                firstNotFullSegment[outSA] += 1;
                firstAvailableOffset[outSA] = 0;
              }                  
            else
              {
                // The following call increments firstNotFullSegment[] and
                // nSeg[], and sets firstAvailableOffset[] = 0
                add_segment(outSA);
              };
        
          };
        
        //        py::array_t<PS,0> & vec = segListPair[outSA][firstNotFullSegment[outSA]];
        //        auto& vec = segListPair[outSA][firstNotFullSegment[outSA]];

        py::buffer_info pseg_info = segListPair[outSA][firstNotFullSegment[outSA]].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto pseg = static_cast<Pstruct<PT>*>(pseg_info.ptr); // Pointer to a structured Numpy array

        //        vec[firstAvailableOffset[outSA]] = item_input;
        auto ip = firstAvailableOffset[outSA];

        // Loop on the tuple to add the new item data
        //      for(item_input::size_type i = 0; i != item_input.size(); i++)
        //      for(py::ssize_t i = 0; i != item_input.size(); i++)
        /*      
                for(unsigned i = 0; i != item_input.size(); i++)
                {
                pseg[ip][i] = item_input[i];
                }
        */
        pseg[ip].set_from_tuple(item_input); // OR use an overloaded .set() function.

        // Compute the full zero-based index of the particle for return
        auto full_index = firstNotFullSegment[outSA]*segmentLength + firstAvailableOffset[outSA];
        //Increment the next available slot for next time
        firstAvailableOffset[outSA] += 1;

        //      int vec = 1;
        //      return py::make_tuple(vec[firstAvailableOffset[outSA]-1], full_index);
        //      return py::make_tuple(pseg[firstAvailableOffset[outSA]-1], full_index);
        // Return the index in this segment, and the full index
        return py::make_tuple(ip, full_index);
      }
    
    //! Get the Numpy array and offset given a full SA index.
    /*!

      Given the full SegmentedArray index of a stored item, return the Numpy array
      (segment) and the offset of the item in the array.  This is used, e.g., to
      provide access to the item in Python.

      \param full_index is the full SA index of the item.
      \return the tuple: (Numpy array, item offset)

    */
    py::tuple get_segment_and_offset(py::ssize_t full_index)
      {

        // (seg, offset) = divmod(i, self.SEGMENTLENGTH)
        auto seg = (py::ssize_t) full_index / segmentLength;
        auto offset = full_index % segmentLength;
      
        auto outSA = outSegmentedArray;
        /* This is for accessing data values in C++
        py::buffer_info pseg_info = segListPair[outSA][seg].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto pseg = static_cast<Pstruct<PT>*>(pseg_info.ptr); // Pointer to a structured Numpy array
        */
        //        return self.segListPair[outSA][seg][offset];
        
        return py::make_tuple(segListPair[outSA][seg], offset);
        
      }

    //! Get a tuple containing the item values at a full SA index
    py::tuple get_as_tuple(py::ssize_t full_index)
      {
        auto seg = (py::ssize_t) full_index / segmentLength;
        auto offset = full_index % segmentLength;
      
        auto outSA = outSegmentedArray;
        py::buffer_info seg_info = segListPair[outSA][seg].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto parray = static_cast<Pstruct<PT>*>(seg_info.ptr); // Pointer to a structured Numpy array
        return parray[offset].as_tuple();
      }

    //! Initialize a loop over the all the array segments.
    /*!

      The loop should call get_next_segment('in') and get_next_out_segment() in the
      iteration loop.

      \param returnDataPtrs: if true, return pointers to the data instead of Numpy arrays
      \return The tuple: (number of items in the first segment of 'in' SA,
                          ref to first segment of 'in' SA,
                          ref to first segment of 'out' SA)
    */
    py::tuple init_inout_loop(bool returnDataPtrs = false)
      {
        // Abbreviations
        auto inSA = inSegmentedArray;
        auto outSA = outSegmentedArray;

        currentSegment[0] = 0;
        currentSegment[1] = 0;

        /* Swap the two SegmentedArrays in the pair, so that the current 'out' array
           becomes the 'in' array.
           The 'in' array appears on the RHS of expressions, i.e., it contains input
           values for an expression.
           The 'out' array appears on the LHS of expressions, i.e., it is the result
           of evaluating an expression. It's like a scratch array: it's initial values
           don't matter, just it's length.
        */

        inSA = (inSA+1) % 2;
        outSA = (outSA+1) % 2;
        
        py::ssize_t lastItem = 0;
        if (firstNotFullSegment[inSA] == 0)
          {
            lastItem = firstAvailableOffset[inSA];
            if (lastItem == 0) return py::make_tuple(0, nullptr, nullptr);  // WILL THIS WORK FOR RETURNING BOTH NUMPY ARRAY AND DATA POINTERS?
          }
        else
          {
            lastItem = segmentLength;
          }

        inSegmentedArray = inSA;
        outSegmentedArray = outSA;

        py::ssize_t segIndex = 0;
        auto inSeg = segListPair[inSA][segIndex];
        auto outSeg = segListPair[outSA][segIndex];

        // Return either the Numpy arrays, or the pointers to the data
        if (returnDataPtrs == true)
          {
            // Pointers to the data:
            // in:
            py::buffer_info inSeg_info = inSeg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto inSegData = static_cast<Pstruct<PT>*>(inSeg_info.ptr); // Pointer to a structured Numpy array
            // out:
            py::buffer_info outSeg_info = outSeg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto outSegData = static_cast<Pstruct<PT>*>(outSeg_info.ptr); // Pointer to a structured Numpy array

            return py::make_tuple(lastItem, inSegData, outSegData);
          }
        else
          {
            // Numpy arrays:
            return py::make_tuple(lastItem, inSeg, outSeg);
          }

        // orig: return py::make_tuple(lastItem, segListPair[inSA][segIndex], segListPair[outSA][segIndex]);
        
      }
    // py::tuple init_inout_loop(bool returnDataPtrs = false) ENDDEF

    /*! Return a tuple with a reference to the next segment of either the 'in' or
        'out' array, and the number of active items.

        \param in_out is a enum, either in_array or out_array.
        \param returnDataPtr: if true, return a pointer to the data instead of a Numpy array
        \return The tuple (number of active items in the next array segment, 
                           reference to next array segment)

    */
    py::tuple get_next_segment(WhichArray in_out, bool returnDataPtr = false)
      {
        unsigned theSA(2); // The value 2 is for initialization; should never occur
                           // below (see assert).
        if (in_out == WhichArray::in_array)
          theSA = inSegmentedArray;
        else
          if (in_out == WhichArray::out_array)
            theSA = outSegmentedArray;
        
        assert(theSA == 0 || theSA == 1);

        currentSegment[theSA] += 1;
        py::ssize_t segIndex = currentSegment[theSA];

        // If the segment index exceeds occupied limit, return nullptr.
        if (segIndex > firstNotFullSegment[theSA]) return py::make_tuple(0, nullptr);
        //  If this is the last segment and it's empty, return nullptr.
        //  ===> This should cause the caller to break out of the loop <===
        //  If it's not empty, return the non-empty items
        py::ssize_t lastItem;
        if (firstNotFullSegment[theSA] == segIndex)
          {
            lastItem = firstAvailableOffset[theSA];
            if (lastItem == 0) return py::make_tuple(0, nullptr);
          }
        else
          {
            lastItem = segmentLength;
          }

        auto seg = segListPair[theSA][segIndex];
        
        // Return either the Numpy array, or the pointer to the data
        
        if (returnDataPtr == true)
          {
            // Pointer to the data:
            py::buffer_info seg_info = seg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto segData = static_cast<Pstruct<PT>*>(seg_info.ptr); // Pointer to a structured Numpy array
            return py::make_tuple(lastItem, segData);
          }
        else
          {
            // Numpy array:
            return py::make_tuple(lastItem, seg);
          }
        //        return py::make_tuple(lastItem, segListPair[theSA][segIndex]);
      }
    // py::tuple get_next_segment(WhichArray in_out, bool returnDataPtr = false) ENDDEF
    
    //! Return a tuple containing a reference to the next segment of the 'out' array.
    /*!

      If this function is called, it assumes that you need space to write on, so it
      will allocate a new segment if we're out of 'out' segments.

      This method looks similar to the push_back() method above, since the 'out'
      array returned is effectively scratch space, i.e., it is only written to.

      A tuple is returned because the return can be either a Numpy array or a data
      pointer. The tuple allows for either one.

      \param returnDataPtr: if true, return a pointer to the data instead of a Numpy array
      \return: A tuple containing a reference to next segment of the 'out' array.

    */
    py::tuple get_next_out_segment(bool returnDataPtr = false)
      //  py::array_t<Pstruct<PT>, 0> get_next_out_segment(bool returnDataPtr = false)
      {
        // Abbreviations
        auto outSA = outSegmentedArray;
        currentSegment[outSA] += 1;
        py::ssize_t segIndex = currentSegment[outSA];

        // If another segment is already available, use it. Otherwise, allocate a new
        // segment.
        if (segIndex < nSeg[outSA])
          {
            firstNotFullSegment[outSA] += 1;
            firstAvailableOffset[outSA] = 0;
          }
        else
          {
            // The following call increments the variables firstNotFullSegment[] and
            // nSeg[], and sets firstAvailableOffset[] = 0
            add_segment(outSA);
          }

        auto seg = segListPair[outSA][segIndex];

        if (returnDataPtr == true)
          {
            // Pointer to the data:
            py::buffer_info seg_info = seg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto segData = static_cast<Pstruct<PT>*>(seg_info.ptr); // Pointer to a structured Numpy array
            return py::make_tuple(segData);
          }
        else
          {
            // Numpy array:
            return py::make_tuple(seg);
          }
        // orig: return segListPair[outSA][segIndex];
      }
    // py::tuple get_next_out_segment(bool returnDataPtr = false) ENDDEF
    
  };
  // class SegmentedArrayPair ENDCLASS


/*
  void pseg_DnT_pstruct_instances()
  {
    SegmentedArrayPair<DnT_pstruct1D> seg_array_obj(1);
  }
*/
  
} // namespace dnt

#endif
