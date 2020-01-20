/*! \file SegmentedArrayPair.h

  \brief This file has the source code for a C++ implementation of SegmentedArrayPair storage.

  SA is an abbreviation for SegmentedArray
  SAP is an abbreviation for SegmentedArrayPair

  \namespace dnt
  \sa segmented_array_pair_cpp.cpp

*/

#ifndef SEGMENTEDARRAYPAIR_H
#define SEGMENTEDARRAYPAIR_H

#include <iostream>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include "Pstruct.h"

namespace dnt
{
  // WhichArray is NOT USED anywhere.
  /*! \class WhichArray

      \brief WhichArray is an enum class to tag the 'I' and 'O' arrays of a
             SegmentedArrayPair

  */
  enum class WhichArray
  {
    //   'I',
    //   'O',
    in_array='I',
    out_array='O'
  };
  
/*! \class SegmentedArrayPair
    \brief The SegmentedArrayPair class is a C++ version of the Python class SegmentedArrayPair_C.

    This class implements a segmented array of items.  Each segment is a Numpy
    structured array, i.e., each item in the array is a struct. The data fields in
    the struct have names and data types. All segments are the same length and use
    the same struct type.

    The class is templated on the struct (PT) that it has to store.

    The segments are always full of structs, except possibly the last one.

    SegList[0] is the first Numpy array.
    SegList[n-1] is the nth Numpy array.

    \param PT is the type of particle struct.
    \sa Ptype, Pstruct

*/
  template<Ptype PT>
    class SegmentedArrayPair
  {

  private:
    // The number of items in one contiguous segment.
    const py::ssize_t segment_length;

    // Make a pair of empty list of segments

    // This contains two lists of Numpy arrays
    std::vector<py::array_t<Pstruct<PT>,0>> seg_list_pair[2];

    unsigned nseg[2];
    unsigned npmax[2];
    
    // Segment number and offset for the first available opening
    unsigned first_not_full_segment[2];
    unsigned first_available_offset[2];
    
    // current_segment is used to count segments in loops over
    // segments, e.g., to push particles.
    unsigned current_segment[2];

    // Identify which of the array pairs is being written to. When
    // the particles are being initialized, we write to the first
    // one of the array pair, i.e., index 0.
    unsigned in_segmented_array;
    unsigned out_segmented_array;


    //! Create a structured Numpy array of length *segment_length*.
    /*!
      Function allocate_a_segment() creates a Numpy array of length *segment_length*
      for the PT particle type, and returns a pointer to it.

      \return pointer to the new Numpy array.

    */
    py::array_t<Pstruct<PT>, 0> allocate_a_segment(void)
    {
      // No pointer is passed, so Numpy will allocate the buffer
      return py::array_t<Pstruct<PT>, 0>(segment_length);
    }
    // ENDDEF: py::array_t<Pstruct<PT>, 0> allocate_a_segment(void)


    //! Add another segment to the specified SegmentedArray in order to store more items.
    /*!

      \param theSA (0 or 1) identifies which of the two SegmentedArrays needs a new segment.
      \return void

    */
    void add_segment(unsigned theSA)
    {
      nseg[theSA] += 1;
      //      py::array_t<Pstruct<PT>, 0> newSegment = allocate_a_segment(segment_length);
      //      py::array_t<Pstruct<PT>, 0> newSegment = allocate_a_segment();
      // No pointer is passed, so Numpy will allocate the buffer
      py::array_t<Pstruct<PT>, 0> newSegment = py::array_t<Pstruct<PT>, 0>(segment_length);
      seg_list_pair[theSA].push_back(newSegment);

      first_not_full_segment[theSA] += 1;
      first_available_offset[theSA] = 0;
    }
    // ENDDEF: void add_segment(unsigned theSA)
    
  public:
    //! The one and only SegmentedArrayPair ctor.
    /*!

      \param segment_length is the number of structs stored in one segment of the array.

      \sa Ptype

    */
    SegmentedArrayPair(py::ssize_t segment_length):
      segment_length(segment_length),
      nseg{0, 0},
      npmax{0, 0},
      first_not_full_segment{0, 0},
      first_available_offset{0, 0},
      current_segment{0, 0},
      in_segmented_array(1),
      out_segmented_array(0)
      {
        //        std::cout << "Hello from the SegmentedArrayPair ctor" << std::endl;
        for (auto iSA : {0, 1})
          {
            // Add the Numpy arrays for the first segment pair
            //            py::array_t<Pstruct<PT>, 0> newSegment = allocate_a_segment(segment_length);
            //            py::array_t<Pstruct<PT>, 0> newSegment = allocate_a_segment();
            py::array_t<Pstruct<PT>, 0> newSegment = py::array_t<Pstruct<PT>, 0>(segment_length);
            seg_list_pair[iSA].push_back(newSegment);
            // Count the number of segments:
            nseg[iSA] = seg_list_pair[iSA].size();
            // Maximum number of particles that can be stored at present
            npmax[iSA] = nseg[iSA]*segment_length;
          }
      }
      // ENDDEF: SegmentedArrayPair(py::ssize_t segment_length)


    // The dtor (See pybind 8.5 Non-public destructors)
    ~SegmentedArrayPair()
      {
        std::cout << "The dtor ~SegmentedArrayPair has been called" << std::endl;
      
        // Release the Numpy arrays?  No: this causes a crash. The
        // syntax may be wrong below. How should this memory be released?
      
        /* for (auto iSA : {0, 1}) */
        /* { */
        /*   for (auto const& seg : seg_list_pair[iSA]) */
        /*   { */
        /*     delete &seg; // Is this correct? */
        /*   }       */
        /* } */
      } 


    //! Return the segment length
    py::ssize_t get_segment_length()
      {
        return segment_length;
      }
    // ENDDEF: py::ssize_t get_segment_length()


    //! Return a py::tuple with the current number of segments in the "in" and "out" arrays.
    py::tuple get_number_of_segments()
      {
        // Abbreviations
        auto inSA = in_segmented_array;
        auto outSA = out_segmented_array;
        
        return py::make_tuple(nseg[inSA], nseg[outSA]);
      }
    // ENDDEF: py::tuple get_number_of_segments()


    // BEGINDEF: py::ssize_t get_number_of_items()
    //! Return the number of items stored in the "out" array
      py::ssize_t get_number_of_items()
      {
        // Abbreviations
        auto outSA = out_segmented_array;

        return first_not_full_segment[outSA]*segment_length + first_available_offset[outSA];
      }
      // ENDDEF: py::ssize_t get_number_of_items()


      //! Sets the number of active items currently stored.
      /*! 
          \param string in_out: Either "in" or "out" depending on whether
                                we're dealing with the in or out SA.

          \return: Nothing is returned.
      */
      void set_number_of_items(std::string in_out, py::ssize_t n_items)
      {
        unsigned theSA(2); // Initialized to avoid a "may be uninitialized" warning
        if (in_out.compare("in") == 0)
          theSA = in_segmented_array;
        else
          if (in_out.compare("out") == 0)
            theSA = out_segmented_array;

        assert(theSA == 0 || theSA == 1);

        // Compute the segment and offset of the last item. The
        // zero-based index of the last item is n_items-1
        // (seg, offset) = divmod(n_items-1, self.segment_length)
        auto seg = (py::ssize_t) (n_items-1) / segment_length;
        auto offset = (n_items - 1) % segment_length;

        first_not_full_segment[theSA] = seg;
        first_available_offset[theSA] = offset+1;

        return;
      }
      // ENDDEF: void set_number_of_items(std::string in_out, py::ssize_t n_items)      


    //! Return the total number of items that can currently be stored.
    py::tuple get_capacity()
      {
        // Abbreviations
        auto inSA = in_segmented_array;
        auto outSA = out_segmented_array;

        return py::make_tuple(nseg[inSA]*segment_length, nseg[outSA]*segment_length);
      }
    // ENDDEF: py::tuple get_capacity()

    
    //! Return the number of megabytes allocated for the item arrays
    py::tuple get_number_of_mbytes()
      {
        // Abbreviations
        auto inSA = in_segmented_array;
        auto outSA = out_segmented_array;

        // Get the number of bytes in one segment (a structured Numpy array)
        py::buffer_info segInfo = seg_list_pair[outSA][0].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        py::ssize_t itembytes = segInfo.itemsize;
        py::ssize_t nitems = segInfo.size;
        py::ssize_t nbytes = nitems*itembytes;
        
        return py::make_tuple(nseg[inSA]*nbytes/(1.0e6), nseg[outSA]*nbytes/(1.0e6));
      }
    // ENDDEF: py::tuple get_number_of_mbytes()    


    //! Add data from a 1-element py::array_t to the "out" SegmentedArray.
    /*!

      \param input_item is a 1-element py::array_t<Pstruct<PT>, 0> containing a complete item.
      \return the tuple: (array idex containing the new item, full SA index)

    */
    py::tuple push_back(py::array_t<Pstruct<PT>, 0> input_item)
      {
        // Abbreviations
        auto outSA = out_segmented_array;

        // Locate a slot for the item.
        /*        
           If we've reached the end of the current segment, we need to switch to the
           next segment, if there is one, or else add a new segment.
        */
        if (first_available_offset[outSA] == segment_length)
          {
            // If another segment is already available, use
            // it. Otherwise, allocate a new segment.
            current_segment[outSA] += 1;
            if (current_segment[outSA] < nseg[outSA])
              {
                first_not_full_segment[outSA] += 1;
                first_available_offset[outSA] = 0;
              }                  
            else
              {
                // The following call increments first_not_full_segment[] and
                // nseg[], and sets first_available_offset[] = 0
                add_segment(outSA);
              };
        
          };

        // Use the buffer protocol to get to the contiguous data in the Numpy arrays
        // "out" SA
        py::buffer_info segInfo = seg_list_pair[outSA][first_not_full_segment[outSA]].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto pseg = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to the contiguous data.
        // See docs above for *input_item*
        py::buffer_info itemInputInfo = input_item.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto inputArray = static_cast<Pstruct<PT>*>(itemInputInfo.ptr); // Pointer to the contiguous data.
        
        // Loop on the tuple to add the new item data
        //      for(input_item::size_type i = 0; i != input_item.size(); i++)
        //      for(py::ssize_t i = 0; i != input_item.size(); i++)
        /*      
                for(unsigned i = 0; i != input_item.size(); i++)
                {
                pseg[ip][i] = input_item[i];
                }
        */
        // Copy the item to the end of the SA
        auto ip = first_available_offset[outSA];
        pseg[ip] = inputArray[0];
        
        // Compute the full zero-based index of the particle for return
        auto full_index = first_not_full_segment[outSA]*segment_length + first_available_offset[outSA];
        //Increment the next available slot for next time
        first_available_offset[outSA] += 1;

        // Return the index into this segment, and the full index
        return py::make_tuple(ip, full_index);
      }
    // ENDDEF: py::tuple push_back(py::array_t<Pstruct<PT>, 0> input_item)


    //! Add data from a py::tuple to the "out" SegmentedArray.
    /*!

      \param item_input is a tuple containing a complete item.
      \return the tuple: (array idex with the new item, full SA index)

    */
    py::tuple push_back(py::tuple item_input)
      {
        // Abbreviations
        auto outSA = out_segmented_array;

        /*        
                  If we've reached the end of the current segment, we need to
                  switch to the next segment, if there is one, or else add a new
                  segment.
        */
        if (first_available_offset[outSA] == segment_length)
          {
            // If another segment is already available, use
            // it. Otherwise, allocate a new segment.
            current_segment[outSA] += 1;
            if (current_segment[outSA] < nseg[outSA])
              {
                first_not_full_segment[outSA] += 1;
                first_available_offset[outSA] = 0;
              }                  
            else
              {
                // The following call increments first_not_full_segment[] and
                // nseg[], and sets first_available_offset[] = 0
                add_segment(outSA);
              };
          };
        
        py::buffer_info segInfo = seg_list_pair[outSA][first_not_full_segment[outSA]].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto pseg = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array

        // Loop on the tuple to add the new item data
        //      for(item_input::size_type i = 0; i != item_input.size(); i++)
        //      for(py::ssize_t i = 0; i != item_input.size(); i++)
        /*      
                for(unsigned i = 0; i != item_input.size(); i++)
                {
                pseg[ip][i] = item_input[i];
                }
        */
        auto ip = first_available_offset[outSA];
        pseg[ip].set_from_list_or_tuple(item_input); // OR use an overloaded .set() function.

        // Compute the full zero-based index of the particle for return
        auto full_index = first_not_full_segment[outSA]*segment_length + first_available_offset[outSA];
        //Increment the next available slot for next time
        first_available_offset[outSA] += 1;

        return py::make_tuple(ip, full_index);
      }
    // ENDDEF: py::tuple push_back(py::tuple item_input)

// #undef LIST
#define LIST
#ifdef LIST 
    //! Add data from a py::list to the "out" SegmentedArray.
    /*!

      \param item_input is a py::list containing a complete item.
      \return the tuple: (array index with the new item, full SA index)

    */
    py::tuple push_back(py::list item_input)
      {
        // Abbreviations
        auto outSA = out_segmented_array;

        /*        
          If we've reached the end of the current segment, we need to switch to the
          next segment, if there is one, or else add a new segment.
        */
        if (first_available_offset[outSA] == segment_length)
          {
            // If another segment is already available, use
            // it. Otherwise, allocate a new segment.
            current_segment[outSA] += 1;
            if (current_segment[outSA] < nseg[outSA])
              {
                first_not_full_segment[outSA] += 1;
                first_available_offset[outSA] = 0;
              }                  
            else
              {
                // The following call increments first_not_full_segment[] and
                // nseg[], and sets first_available_offset[] = 0
                add_segment(outSA);
              };
        
          };
        
        py::buffer_info segInfo = seg_list_pair[outSA][first_not_full_segment[outSA]].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto pseg = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array

        auto ip = first_available_offset[outSA];
        pseg[ip].set_from_list_or_tuple(item_input); // OR use an overloaded .set() function.

        // Compute the full zero-based index of the particle for return
        auto full_index = first_not_full_segment[outSA]*segment_length + first_available_offset[outSA];
        //Increment the next available slot for next time
        first_available_offset[outSA] += 1;

        // Return the index in this segment, and the full index
        return py::make_tuple(ip, full_index);
      }
    // ENDDEF: py::tuple push_back(py::list item_input)
#endif

    
    //! Get the Numpy array and offset given a full SA index.
    /*!

      Given the full SegmentedArray index of a stored item, return a reference to the
      Numpy array (segment) and the offset of the item in the array.  This can be
      used, e.g., to provide access to the stored item in Python.

      \param full_index is the full SA index of the item.
      \return the tuple: (Numpy array reference, item offset)

    */
    py::tuple get_segment_and_offset(py::ssize_t full_index)
      {
        auto seg = (py::ssize_t) full_index / segment_length;
        auto offset = full_index % segment_length;
      
        auto outSA = out_segmented_array;
        
        return py::make_tuple(seg_list_pair[outSA][seg], offset);
      }
    // ENDDEF: py::tuple get_segment_and_offset(py::ssize_t full_index)


    //! Return a copy of an item, as a py::dict, given the full SA index
    py::dict get_item(py::ssize_t full_index)
      {
        auto seg = (py::ssize_t) full_index / segment_length;
        auto offset = full_index % segment_length;
      
        auto outSA = out_segmented_array;
        py::buffer_info segInfo = seg_list_pair[outSA][seg].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto parray = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array
        // parray[offset]; // contains the particle data
          // Make a dictionary from this
        return parray[offset].as_dict();  // as_dict() is defined in Pstruct.h
      }
    // ENDDEF: py::dict get_item(py::ssize_t full_index)
    
    
    //! Get a copy of an item's values, as a py::tuple, given the full SA index of the item.
    /*!  Note that the returned py::tuple cannot be modified. If the returned copy of
         the item values is to be modified, use get_as_list() instead.

         If the stored item needs to be modified, use get_segment_and_offset() to get
         access to it.

         \sa get_segment_and_offset, get_as_list
     */
    py::tuple get_as_tuple(py::ssize_t full_index)
      {
        auto seg = (py::ssize_t) full_index / segment_length;
        auto offset = full_index % segment_length;
      
        auto outSA = out_segmented_array;
        py::buffer_info segInfo = seg_list_pair[outSA][seg].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto parray = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array
        return parray[offset].as_tuple();
      }
    // ENDDEF: py::tuple get_as_tuple(py::ssize_t full_index)


    //! Get a copy of an item's values, as a py::list, given the full SA index of the item.
    /*!  The returned py::list can be modified, but this does not modify the stored item.

      If the stored item needs to be modified, use get_segment_and_offset() to get
      access to it.

      \sa get_segment_and_offset, get_as_tuple

    */
    py::list get_as_list(py::ssize_t full_index)
      {
        auto seg = (py::ssize_t) full_index / segment_length;
        auto offset = full_index % segment_length;
      
        auto outSA = out_segmented_array;
        py::buffer_info segInfo = seg_list_pair[outSA][seg].request(); // request() returns metadata about the array (ptr, ndim, size, shape)
        const auto parray = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array
        return parray[offset].as_list();
      }
    // ENDDEF: py::list get_as_list(py::ssize_t full_index)


    //! Initialize a loop over the segments of the "out" array.
    /*!

      This function is used to start a loop that doesn't change the positions of the
      items, so the "in" array isn't needed.

      \param return_data_ptr: if true, return pointers to the data instead of Numpy arrays
      \return The tuple: (number of items in the first segment of "out" SA,
                          ref to first segment of "out" SA)
    */
    py::tuple init_out_loop(bool return_data_ptr = false)
      //    py::tuple init_out_loop()
      {
        // Abbreviations
        auto outSA = out_segmented_array;

        // std::cout << "init_out_loop(): out segment is using pair index " << outSA << std::endl;
        
        // This is used to count through the segments
        current_segment[outSA] = 0;

        py::ssize_t lastItem = 0;
        // If the first segment is also the last segment, the number of items in it
        // may be less than the length of the segment
        if (first_not_full_segment[outSA] == 0)
          {
            lastItem = first_available_offset[outSA];
            if (lastItem == 0) return py::make_tuple(0, nullptr);  // WILL THIS WORK FOR RETURNING BOTH NUMPY ARRAY AND DATA POINTERS?
          }
        else
          {
            lastItem = segment_length;
          }

        py::ssize_t segIndex = 0;
        auto outSeg = seg_list_pair[outSA][segIndex];

        // Return either the Numpy arrays, or the pointers to the data
        // bool return_data_ptr = false;
        if (return_data_ptr == true)
          {
            // Pointers to the data:
            py::buffer_info outSegInfo = outSeg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto outSegData = static_cast<Pstruct<PT>*>(outSegInfo.ptr); // Pointer to a structured Numpy array
            return py::make_tuple(lastItem, outSegData);
          }
        else
          {
            // Numpy arrays:
            return py::make_tuple(lastItem, outSeg);
          }
      }
    // ENDDEF: py::tuple init_out_loop(bool return_data_ptr = false)


    bool swapPair = true; // Set to false if the 'in' and 'out' SAs should not be
                          // swapped. For testing, uncomment the swapPair =! swapPair
                          // statement below. Don't modify this statement.
    //! Initialize a loop over the all the array segments.
    /*!

      The loop should call get_next_segment("in") and get_next_out_segment() in the
      iteration loop.

      \param return_data_ptrs: if true, return pointers to the data instead of Numpy arrays
      \return The tuple: (number of items in the first segment of "in" SA,
      ref to first segment of "in" SA,
      ref to first segment of "out" SA)
    */
    py::tuple init_inout_loop(bool return_data_ptrs = false)
      {
        // Abbreviations
        auto inSA = in_segmented_array;
        auto outSA = out_segmented_array;

        current_segment[0] = 0;
        current_segment[1] = 0;

        /* Swap the two SegmentedArrays in the pair, so that the current "out" array
           becomes the "in" array.
           The "in" array appears on the RHS of expressions, i.e., it contains input
           values for an expression.
           The "out" array appears on the LHS of expressions, i.e., it is the result
           of evaluating an expression. It's like a scratch array: it's initial values
           don't matter, just it's length.
        */
        
        if (swapPair)
          {
            inSA = (inSA+1) % 2;
            outSA = (outSA+1) % 2;
          }

        /*
          Uncomment the following for testing only. E.g., suppose we call this function,
          but the subsequent code don't actually copy the "in" to the "out" segment, and
          then we want to call this function again. By setting swapPair = false, we will
          get back the same segments as on the first call. Otherwise, the segments will be
          swapped on the second call.
        */
        // swapPair = !swapPair;
        
        py::ssize_t lastItem = 0;
        if (first_not_full_segment[inSA] == 0)
          {
            lastItem = first_available_offset[inSA];
            if (lastItem == 0) return py::make_tuple(0, nullptr, nullptr);  // WILL THIS WORK FOR RETURNING BOTH NUMPY ARRAY AND DATA POINTERS?
          }
        else
          {
            lastItem = segment_length;
          }

        in_segmented_array = inSA;
        out_segmented_array = outSA;

        py::ssize_t segIndex = 0;
        auto inSeg = seg_list_pair[inSA][segIndex];
        auto outSeg = seg_list_pair[outSA][segIndex];

        // std::cout << "init_inout_loop(): out segment is using pair index " << outSA << std::endl;

        // Return either the Numpy arrays, or the pointers to the data
        // bool return_data_ptrs = false;
        if (return_data_ptrs == true)
          {
            // Pointers to the data:
            // "in" array:
            py::buffer_info inSegInfo = inSeg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto inSegData = static_cast<Pstruct<PT>*>(inSegInfo.ptr); // Pointer to a structured Numpy array
            // "out" array:
            py::buffer_info outSegInfo = outSeg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto outSegData = static_cast<Pstruct<PT>*>(outSegInfo.ptr); // Pointer to a structured Numpy array

            return py::make_tuple(lastItem, py::capsule(inSegData, "inSegData pointer", nullptr), py::capsule(outSegData, "outSegData pointer", nullptr));
          }
        else
          {
            // Numpy arrays:
            return py::make_tuple(lastItem, inSeg, outSeg);
          }
      }
    // ENDDEF: py::tuple init_inout_loop(bool return_data_ptrs = false)


    //! Return a tuple with a reference to the next segment of either the "in" or "out" array, and the number of active items.
    /*      
        \param in_out is a string, either "in" or "out".
        \param return_data_ptr: if true, return a pointer to the data instead of a Numpy array
        \return The tuple (number of active items in the next array segment, 
                           reference to next array segment)

    */
    py::tuple get_next_segment(std::string in_out, bool return_data_ptr = false)
      {
        unsigned theSA(2); // The value 2 is for initialization; should never occur
                           // below (see assert).
        if (in_out.compare("in") == 0)
          theSA = in_segmented_array;
        else
          if (in_out.compare("out") == 0)
            theSA = out_segmented_array;
        
        assert(theSA == 0 || theSA == 1);

        current_segment[theSA] += 1;
        py::ssize_t segIndex = current_segment[theSA];

        // If the segment index exceeds occupied limit, return nullptr.
        if (segIndex > first_not_full_segment[theSA]) return py::make_tuple(0, nullptr);

        //  If this is the last segment and it's empty, return nullptr.
        //  ===> This should cause the caller to break out of the loop <===
        //  If it's not empty, return the non-empty items
        py::ssize_t lastItem;
        if (first_not_full_segment[theSA] == segIndex)
          {
            lastItem = first_available_offset[theSA];
            if (lastItem == 0) return py::make_tuple(0, nullptr);
          }
        else
          {
            lastItem = segment_length;
          }

        auto seg = seg_list_pair[theSA][segIndex];
        
        // Return either the Numpy array, or the pointer to the data
        // bool return_data_ptr = false;        
        if (return_data_ptr == true)
          {
            // Pointer to the data:
            py::buffer_info segInfo = seg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto segData = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array
            //            return py::make_tuple(lastItem, segData);
            return py::make_tuple(lastItem, py::capsule(segData));
          }
        else
          {
            // Numpy array:
            return py::make_tuple(lastItem, seg);
          }
        //        return py::make_tuple(lastItem, seg_list_pair[theSA][segIndex]);
      }
    // ENDDEF py::tuple get_next_segment(std::string in_out, bool return_data_ptr = false)

    
    //! Return a tuple containing a reference to the next segment of the "out" array.
    /*!

      When this function is called, it assumes that you need space to write on, so it
      will allocate a new segment if we're out of "out" segments.

      This method looks similar to the push_back() method above, since the "out"
      array returned is effectively scratch space, i.e., it is only written to.

      A tuple is returned because the return can be either a Numpy array or a data
      pointer. The tuple allows for either one.

      \param return_data_ptr: if true, return a pointer to the data instead of a Numpy array
      \return: A tuple containing a reference to next segment of the "out" array.

    */
    py::tuple get_next_out_segment(bool return_data_ptr = false)
      {
        // Abbreviations
        auto outSA = out_segmented_array;
        current_segment[outSA] += 1;
        py::ssize_t segIndex = current_segment[outSA];

        // If another segment is already available, use it. Otherwise, allocate a new
        // segment.
        if (segIndex < nseg[outSA])
          {
            first_not_full_segment[outSA] += 1;
            first_available_offset[outSA] = 0;
          }
        else
          {
            // The following call increments the variables first_not_full_segment[] and
            // nseg[], and sets first_available_offset[] = 0
            add_segment(outSA);
          }

        auto seg = seg_list_pair[outSA][segIndex];
        // bool return_data_ptr = false;
        if (return_data_ptr == true)
          {
            // Pointer to the data:
            py::buffer_info segInfo = seg.request(); // request() returns metadata about the array (ptr, ndim, size, shape)
            const auto segData = static_cast<Pstruct<PT>*>(segInfo.ptr); // Pointer to a structured Numpy array
            return py::make_tuple(segData);
          }
        else
          {
            // Numpy array:
            return py::make_tuple(seg);
          }
      }
    // ENDDEF: py::tuple get_next_out_segment(bool return_data_ptr = false)

    
    //! Return a reference to the current segment of the "out" array.
    /*!

      \return: A reference to current segment of the "out" array.

    */
    py::array_t<Pstruct<PT>,0> get_current_out_segment(void)
      {
        // Abbreviations
        auto outSA = out_segmented_array;
        py::ssize_t segIndex = current_segment[outSA];

        return seg_list_pair[outSA][segIndex];
      }
    // ENDDEF: py::tuple get_current_out_segment(bool return_data_ptr = false)
                      
    //! Get the full index of an item in either the "in" or "out" SA.
    /*!  
      \param indx[in] offset of an item into the SA.
      \param in_out[in] "in" or "out", depending of which SA is intended.
      \return the full index of the item.
    */
    py::ssize_t get_full_index(py::ssize_t indx, std::string in_out)
      {
        unsigned theSA(2); // Initialized to avoid a "may be uninitialized" warning
        if (in_out.compare("in") == 0)
          theSA = in_segmented_array;
        else
          if (in_out.compare("out") == 0)
            theSA = out_segmented_array;

        assert(theSA == 0 || theSA == 1);

        py::ssize_t full_index = current_segment[theSA]*segment_length + indx;
        return full_index;
      }
    // ENDDEF: py::ssize_t get_full_index(py::ssize_t indx, std::string in_out)
    
                      
    //! Get the full indices of an item in the "in" and "out" SA.
    /*!  
      \param i_in: offset of an item into the "in" SA.
      \param i_out: offset of an item into the "out" SA.

      \return the tuple (full index in the "in" array, full index in the "out" array)
    */
    py::tuple get_full_indices(py::ssize_t i_in, py::ssize_t i_out)
      {
        // Abbreviations
        auto inSA = in_segmented_array;
        auto outSA = out_segmented_array;

        py::ssize_t full_index_in = current_segment[inSA]*segment_length + i_in;
        py::ssize_t full_index_out = current_segment[outSA]*segment_length + i_out;

        return py::make_tuple(full_index_in, full_index_out);
      }
    // ENDDEF: py::tuple get_full_indices(py::ssize_t i_in, py::ssize_t i_out)

  };
  // class SegmentedArrayPair ENDCLASS

} // namespace dnt

#endif
