/*! \file SegmentedArray_Cpp.h

  \brief This file has the source code for a C++ implementation of SegmentedArray storage.

  \namespace dnt

  \sa SegmentedArrayPair_Cpp
*/

#ifndef SEGMENTEDARRAYPAIR_CPP_H
#define SEGMENTEDARRAYPAIR_CPP_H

// C++ program to demonstrate factory method design pattern 
#include <iostream>
#include "pstruct.h"

//using namespace std; 

//template<typename PS>
//using SegList = std::vector<py::array_t<PS,0>>;

namespace dnt
{
/*! \class SegmentedArrayPair_Cpp
    \brief The SegmentedArrayPair_Cpp class is a C++ version of the Python class SegmentedArrayPair_C.

    This class implements a segmented array of items.  Each segment is
    a Numpy structured array, i.e., each item in the array is a
    struct. The data fields in the struct have names and data
    types. All segments are the same length and use the same struct
    type.

    The class is templated on the struct (PS) that it has to
    store.

    The segments are always full of structs, except possibly the
    last one.

    SegList[0] is the first Numpy array.
    SegList[n-1] is the nth Numpy array.

    \sa Ptype, DnT_pstruct1D, DnT_pstruct2D, DnT_pstruct3D

*/
//  template<typename T>
  template<Ptype PT>
  class SegmentedArrayPair_Cpp
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

      \param segmentLength is the array length
      \return pointer to the Numpy array

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

    //    py::ssize_t segmentLength;
    //    std::vector<py::array_t<Pstruct<PT>,0>> segListPair[2];    
    //    unsigned outSegmentedArray;
    
    //! The one and only ctor.
    /*!

      \param PT is the type of particle struct.
      \param segment_length is the number of structs stored in one segment of the array.

      \sa Ptype

    */
  SegmentedArrayPair_Cpp(py::ssize_t segment_length): segmentLength(segment_length), nSeg{0, 0}, nPmax{0, 0}, firstNotFullSegment{0, 0}, firstAvailableOffset{0, 0}, currentSegment{0, 0}, inSegmentedArray(1), outSegmentedArray(0) {
      std::cout << "Hello from the SegmentedArrayPair_Cpp ctor" << std::endl;
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
    ~SegmentedArrayPair_Cpp()
      {
      
        std::cout << "The dtor ~SegmentedArrayPair_Cpp has been called" << std::endl;
      
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
    
    //! Return a py::tuple with the number of segments in the two segment lists.
    py::tuple get_number_of_segments()
      {
        return py::make_tuple(nSeg[0], nSeg[1]);
      }

    py::ssize_t get_number_of_segments0()
      //    int get_number_of_segments()
      {
        return nSeg[0];
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
      and the offset of the item in the array.  This is used, e.g., to provide access
      to the item in Python.

      \param full_index is the full SA index of the item.
      \return the tuple: (Numpy array, item offset)

    */
    py::tuple get_array_and_offset(py::ssize_t full_index)
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
    
  }; // class SegmentedArrayPair_Cpp


/*
  void pseg_DnT_pstruct_instances()
  {
    SegmentedArrayPair_Cpp<DnT_pstruct1D> seg_array_obj(1);
  }
*/
  
} // namespace dnt

#endif
