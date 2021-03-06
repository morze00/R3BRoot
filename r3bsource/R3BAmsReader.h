#ifndef R3BAMSREADER_H
#define R3BAMSREADER_H

#include "R3BReader.h"

class TClonesArray;

struct EXT_STR_h101_AMS_t;
typedef struct EXT_STR_h101_AMS_t EXT_STR_h101_AMS;
class ext_data_struct_info;

/**
 * A reader of AMS data with UCESB.
 * Receives mapped raw data and converts it to R3BRoot objects.
 * @author J.L. Rodriguez
 * @since May 12, 2018
 */
class R3BAmsReader : public R3BReader {
	public:
		/**
		 * Standard constructor.
		 * Creates instance of the reader. To be called in the steering macro.
		 * @param Pointer to a full C structure generated by the Ucesb unpacker.
		 */
		R3BAmsReader(EXT_STR_h101_AMS *, UInt_t);
		~R3BAmsReader();

		Bool_t Init(ext_data_struct_info *);
		Bool_t Read();
		void Reset();

                /** Accessor to select online mode **/
                void SetOnline(Bool_t option){fOnline=option;} 

	private:
		/* An event counter */
		unsigned int fNEvent;
		/* Reader specific data structure from ucesb */
		EXT_STR_h101_AMS* fData;
		/* Data offset */
		UInt_t fOffset;
                //Don't store data for online
                Bool_t fOnline;
                /**< Output array. */
	        TClonesArray* fArray; 

	public:
		ClassDef(R3BAmsReader, 0);
};

#endif
