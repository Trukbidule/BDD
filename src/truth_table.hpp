#pragma once

//make_signal( index_t index, bool complement = false )
//get_index( signal_t signal )
//set_complemented(signal_t signal, bool complement)
//is_complemented( signal_t signal )
//get_node( signal_t signal )


#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <cmath>

//nb of bit per block in vector bits
#define BLOCK_SIZE 64

/* masks used to filter out unused bits */
static const uint64_t length_mask[] = {
  0x0000000000000001,
  0x0000000000000003,
  0x000000000000000f,
  0x00000000000000ff,
  0x000000000000ffff,
  0x00000000ffffffff,
  0xffffffffffffffff};

/* masks used to get the bits where a certain variable is 1 */
static const uint64_t var_mask_pos[] = {
  0xaaaaaaaaaaaaaaaa,
  0xcccccccccccccccc,
  0xf0f0f0f0f0f0f0f0,
  0xff00ff00ff00ff00,
  0xffff0000ffff0000,
  0xffffffff00000000};

/* masks used to get the bits where a certain variable is 0 */
static const uint64_t var_mask_neg[] = {
  0x5555555555555555,
  0x3333333333333333,
  0x0f0f0f0f0f0f0f0f,
  0x00ff00ff00ff00ff,
  0x0000ffff0000ffff,
  0x00000000ffffffff};

//return 2^n
uint32_t power_two(const uint32_t n ){
    uint32_t val = 1;
    for(int i=0; i<n; i++){
        val *=2;
    }
    return val;
}

/* return i if n == 2^i*/
inline uint8_t inv_power_two( const uint32_t n ){
    float value = log2(n);
    uint8_t ret= 0u;
    
    //if integer, return the value
    if(ceilf(value) == value)
        ret = (uint8_t)(ceilf(value));
    
    return ret;
}

class Truth_Table{
public:
    //default constructor with thruth table = polarity (all 0, or 1)
    Truth_Table(uint8_t p_num_var) : num_var(p_num_var){
        //this->num_var = p_num_var;
        uint32_t nb_block = (power_two(p_num_var)/BLOCK_SIZE) +1;
        //build the blocks
        for(int i=0; i<nb_block; i++){
            this->bits.push_back(0u);
        }
    }

    //constructor from nb var and truth table in vector
    Truth_Table(uint8_t p_num_var, std::vector<uint64_t> p_bits) : num_var(p_num_var){
        //this->num_var = p_num_var;
        uint32_t nb_block = (power_two(p_num_var)/BLOCK_SIZE) +1;
        uint64_t bit_buffer = 0u;
        
        //safety check: must give enough bits
        assert(nb_block <= p_bits.size());

        for(int i=0; i<nb_block; i++){
            //get the corresponding bits and screen out unecessary bits
            bit_buffer = p_bits.at(i);
            bit_buffer &= mask_of_block(i);
            this->bits.push_back(bit_buffer);//add the bits at index/block i
        }
        
    }

    //constructor from a string of 0 and 1
    Truth_Table(const std::string str) : num_var( inv_power_two(str.size()) ){
        //this->num_var = inv_power_two(str.size());
        uint32_t nb_block = (str.size()/BLOCK_SIZE) +1;

        //safety check
        if(this->num_var == 0)
            return;

        //init the blocks
        for(int i=0; i<nb_block; i++){
            this->bits.push_back(0u);
        }
        
        //check each char, add 1 if find it
        for ( auto i = 0u; i < str.size(); i++ ){
            if ( str[i] == '1' ){
                set_bit( str.size() - 1 - i );
            }else{
                assert( str[i] == '0' );
            }
        }
    }

    //return the value of the bit at the global pos position
    bool get_bit(uint32_t const position) const{
        uint32_t block = position/BLOCK_SIZE;
        uint32_t local_pos = position % BLOCK_SIZE;
        //uint32_t local_pos = position - block*BLOCK_SIZE;
        //assert( position < ( 1 << num_var ) );
        //in the corresponding block, shift the bit at pos, read it
        return ( ( this->bits.at(block) >> local_pos ) & 0x1 );
    }

    void set_bit(uint32_t const position){
        uint32_t block = position/BLOCK_SIZE;
        uint32_t local_pos = position % BLOCK_SIZE;
        
        //std::cout << "block: "<<block <<" local_pos: "<<local_pos<<std::endl;
        //std::cout << "size"<<this->bits.size()<<std::endl;
        //assert( position < ( 1 << num_var ) );
        this->bits.at(block) |= (uint64_t(1) << local_pos);//set the bit
        this->bits.at(block) &= this->mask_of_block(block);//safety: mask unused bits
    }

    //getter num_var
    uint8_t n_var() const{
        return num_var;
    }

    //return the length_mask of the block block
    const uint64_t mask_of_block(const uint32_t block){
        uint32_t last_block = (power_two(this->num_var)/BLOCK_SIZE);
        if(block != last_block)//if not last block: all bit selected
            return length_mask[6];

        //if last block: get the nb of bit in local block
        uint32_t nb_bit_local = power_two(this->num_var) - last_block*BLOCK_SIZE;

        //build a mask by adding 1 shifted in the right positions
        uint32_t mask = 0x0;
        for(int i=0; i<nb_bit_local; i++){
            mask |= (0x1 << i);
        }
        return mask;
    }

    //defined as inline later, after the operators, not asked here
    // Truth_Table positive_cofactor( uint8_t const var ) const;
    // Truth_Table negative_cofactor( uint8_t const var ) const;
    // Truth_Table derivative( uint8_t const var ) const;
    // Truth_Table consensus( uint8_t const var ) const;
    // Truth_Table smoothing( uint8_t const var ) const;

public:
    uint8_t const num_var; /* number of variables involved in the function */
    //uint64_t bits; /* old truth table */
    std::vector<uint64_t> bits; /* the truth table */
};

/*  Reminder for later:
    inline: when fct is stored in memory, copy the inline fct
    within it, to avoid the delay of searching the memory block of the 
    fct called (copied like a macro, not as new fct memorywise)
*/

//block wise NOT: perform NOT operation on multiple blocks in vectors
inline std::vector<uint64_t> bw_NOT(std::vector<uint64_t> v){
    std::vector<uint64_t> v_not;
    
    //NOT on each block, put it in v_not
    for(int i=0; i<v.size(); i++){
        v_not.emplace_back(~v.at(i));
    }
    return v_not;
}

//block wise AND: perform AND operation on multiple blocks in vectors
inline std::vector<uint64_t> bw_AND(std::vector<uint64_t> v1, std::vector<uint64_t> v2){
    assert(v1.size() == v2.size());//safety check
    std::vector<uint64_t> v_and;
    
    //AND on each block, put it in v_and
    for(int i=0; i<v1.size(); i++){
        v_and.emplace_back(v1.at(i) & v2.at(i));
    }
    return v_and;
}

//block wise OR: perform OR operation on multiple blocks in vectors
inline std::vector<uint64_t> bw_OR(std::vector<uint64_t> v1, std::vector<uint64_t> v2){
    assert(v1.size() == v2.size());//safety check
    std::vector<uint64_t> v_or;
    
    //OR on each block, put it in v_or
    for(int i=0; i<v1.size(); i++){
        v_or.emplace_back(v1.at(i) | v2.at(i));
    }
    return v_or;
}

//block wise XOR: perform XOR operation on multiple blocks in vectors
inline std::vector<uint64_t> bw_XOR(std::vector<uint64_t> v1, std::vector<uint64_t> v2){
    assert(v1.size() == v2.size());//safety check
    std::vector<uint64_t> v_xor;
    
    //XOR on each block, put it in v_xor
    for(int i=0; i<v1.size(); i++){
        v_xor.emplace_back(v1.at(i) ^ v2.at(i));
    }
    return v_xor;
}

//block wise equality check
inline bool bw_EQ(std::vector<uint64_t> v1, std::vector<uint64_t> v2){
    if(v1.size() != v2.size())//check block size
        return false;
    
    //check each block
    for(int i=0; i<v1.size(); i++){
        if(v1.at(i) != v2.at(i) )
            return false;
    }
    return true;//if reach the end: everything is equals
}

/* overload std::ostream operator for convenient printing */
inline std::ostream& operator<<(std::ostream& os, Truth_Table const& tt){
    for (int32_t i = (1 << tt.num_var) - 1; i >= 0; --i){
        os << (tt.get_bit(i) ? '1' : '0');
    }
    return os;
}

/* bit-wise NOT operation */
inline Truth_Table operator~(Truth_Table const& tt){
    return Truth_Table(tt.num_var, bw_NOT(tt.bits) );
}

/* bit-wise OR operation */
inline Truth_Table operator|(Truth_Table const& tt1, Truth_Table const& tt2){
    assert( tt1.num_var == tt2.num_var );
    return Truth_Table(tt1.num_var, bw_OR(tt1.bits, tt2.bits) );
}

/* bit-wise AND operation */
inline Truth_Table operator&(Truth_Table const& tt1, Truth_Table const& tt2){
    assert( tt1.num_var == tt2.num_var );
    return Truth_Table( tt1.num_var, bw_AND(tt1.bits, tt2.bits));
}

/* bit-wise XOR operation */
inline Truth_Table operator^(Truth_Table const& tt1, Truth_Table const& tt2){
  assert( tt1.num_var == tt2.num_var );
  return Truth_Table( tt1.num_var, bw_XOR(tt1.bits, tt2.bits) );
}

/* check if two truth_tables are the same */
inline bool operator==(Truth_Table const& tt1, Truth_Table const& tt2){
    if ( tt1.num_var != tt2.num_var ){
        return false;
    }
    return bw_EQ(tt1.bits, tt2.bits);
}

inline bool operator!=(Truth_Table const& tt1, Truth_Table const& tt2){
    if ( tt1.num_var != tt2.num_var ){
        return true;
    }
    return !bw_EQ(tt1.bits, tt2.bits);
}

// inline Truth_Table Truth_Table::positive_cofactor( uint8_t const var ) const{
//     assert( var < num_var );
//     return Truth_Table( num_var, ( bits & var_mask_pos[var] ) | ( ( bits & var_mask_pos[var] ) >> ( 1 << var ) ) );
// }
// 
// inline Truth_Table Truth_Table::negative_cofactor( uint8_t const var ) const{
//     assert( var < num_var );
//     return Truth_Table( num_var, ( bits & var_mask_neg[var] ) | ( ( bits & var_mask_neg[var] ) << ( 1 << var ) ) );
// }
// 
// inline Truth_Table Truth_Table::derivative( uint8_t const var ) const{
//     assert( var < num_var );
//     return positive_cofactor( var ) ^ negative_cofactor( var );
// }
// 
// inline Truth_Table Truth_Table::consensus( uint8_t const var ) const{
//     assert( var < num_var );
//     return positive_cofactor( var ) & negative_cofactor( var );
// }
// 
// inline Truth_Table Truth_Table::smoothing( uint8_t const var ) const{
//     assert( var < num_var );
//     return positive_cofactor( var ) | negative_cofactor( var );
// }

/* Returns the truth table of f(x_0, ..., x_num_var) = x_var (or its complement). */
//generate sequence on the fly, could be optimized to work block by block
Truth_Table create_tt_nth_var( uint8_t const num_var, uint8_t const var, bool const polarity = true ){
    assert (var < num_var );
    
    uint32_t nb_tot_bit = power_two(num_var);
    uint32_t sequence_size = power_two(var);
    
    uint32_t pos_in_sequence = 0;
    bool pol = !polarity;//starting val in the array
    
    Truth_Table tt = Truth_Table(num_var);//empty truth table
    //fill with groups of alternating 0/1 with a period of pow2(var) the tt
    for(uint32_t i=0; i<nb_tot_bit; i++){
        if(pos_in_sequence >= sequence_size){//if starts new sequence
            pos_in_sequence=0;
            pol = !pol;
        }
        if(pol){
            tt.set_bit(i);
        }
        
        pos_in_sequence++;
    }
    
    return tt;
    //return Truth_Table( num_var, polarity ? var_mask_pos[var] : var_mask_neg[var] );
}
