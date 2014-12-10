//--Author	JRM Annand   1st Feb 2007
//--Rev 	JRM Annand
//--Rev 	JRM Annand  11th Jan 2008  1st working version
//--Rev 	JRM Annand...3rd Jun 2008..const Char_t*...gcc 4.3
//--Rev 	JRM Annand..28th Apr 2009..remove TDAQmemmap.h
//--Rev 	JRM Annand..24th Aug 2012..Start debugging
//--Rev 	JRM Annand..25th Sep 2012.."Working version"
//--Update	JRM Annand..25th Sep 2013..Edge detection control
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_V1190
// CAEN VMEbus 32-channel QDC
// Also forms basis for classes to handle similar CAEN ADC hardware
//

#include "TVME_V1190.h"
#include "TDAQexperiment.h"
//
// Internal register offsets
VMEreg_t V1190reg[] = {
  {0x1000,  0x0, 'w', 0x0},          // Control register
  {0x1002,  0x0, 'w', 0x0},          // Status register
  {0x100a,  0x0, 'w', 0x0},          // IRQ level
  {0x100c,  0x0, 'w', 0x0},          // IRQ vector
  {0x100e,  0x0, 'w', 0x0},          // Geo. address reg.
  {0x1010,  0x0, 'w', 0x0},          // MCST base address
  {0x1012,  0x0, 'w', 0x0},          // MCST Control
  {0x1014,  0x0, 'w', 0x0},          // Module reset
  {0x1016,  0x0, 'w', 0x0},          // Software clear
  {0x1018,  0x0, 'w', 0x0},          // Software event reset
  {0x101a,  0x0, 'w', 0x0},          // Software trigger
  {0x101c,  0x0, 'l', 0x0},          // Event counter
  {0x1020,  0x0, 'w', 0x0},          // Event stored
  {0x1022,  0x0, 'w', 0x0},          // Almost full level
  {0x1024,  0x0, 'w', 0x0},          // BLT event number
  {0x1026,  0x0, 'w', 0x0},          // Firmware revision
  {0x1028,  0x0, 'l', 0x0},          // Test register
  {0x102c,  0x0, 'w', 0x0},          // Output control
  {0x102e,  0x0, 'w', 0x0},          // Micro controller data
  {0x1030,  0x0, 'w', 0x0},          // Micro controller handshake
  {0x1032,  0x0, 'w', 0x0},          // Select flash
  {0x1034,  0x0, 'w', 0x0},          // Flash memory
  {0x1036,  0x0, 'w', 0x0},          // SRAM page
  {0x1038,  0x0, 'l', 0x0},          // Event FIFO
  {0x103c,  0x0, 'w', 0x0},          // Event FIFO stored
  {0x103e,  0x0, 'w', 0x0},          // Event FIFO status
  {0x1200,  0x0, 'l', 0x0},          // Dummy 32
  {0x1204,  0x0, 'w', 0x0},          // Dummy 16
  {0x4034,  0x0, 'w', 0x0},          // ID0
  {0x4038,  0x0, 'w', 0x0},          // ID1
  {0x403c,  0x0, 'w', 0x0},          // ID2
  {0x8000,  0x0, 'w', 0x0},          // Compensation SRAM
  {0x0000,  0x0, 'l', 1023},         // output buffer
  {0xffffffff,0x0,'w',0x0},          // terminate
};


enum { EV1190_Threshold=200, EV1190_EdgeDet, EV1190_Window};
static Map_t k1190Keys[] = {
  {"Threshold:",          EV1190_Threshold},
  {"Edge-Detect:",        EV1190_EdgeDet},
  {"Window:",             EV1190_Window},
  {NULL,                  -1}
};

//-----------------------------------------------------------------------------
TVME_V1190::TVME_V1190( Char_t* name, Char_t* file, FILE* log,
			Char_t* line ):
  TVMEmodule( name, file, log, line )
{
  // SetDetection
  // 0 pair, 1 trailing edge, 2 leading edge, 3 leading and trailing
  // Basic initialisation 
  fCtrl = NULL;                            // no control functions
  fType = EDAQ_ADC;                        // analogue to digital converter
  fNreg = EV1190_OutBuff + 1026; // Last "hard-wired" register
  fHardID = 0xa60400;                      // ID read from hardware
  fNBits = 15;                             // 16-bit ADC
  AddCmdList( k1190Keys );                 // CAEN-specific setup commands
  fEdgeDet = 0x2;                          // leading edge
  fWindow = 0x28;
  fWindowOffset = 0xffec;
  fIsContStore = kFALSE;                   // continuous storage
}

//-----------------------------------------------------------------------------
TVME_V1190::~TVME_V1190( )
{
  // Clean up
}

//-----------------------------------------------------------------------------
void TVME_V1190::SetConfig( Char_t* line, Int_t key )
{
  // Configuration from file 
  switch( key ){
  case EV1190_EdgeDet:
    // login input edge detection mode
    // 0 = pair, 1 trailing edge, 2 leading edge, 3 trailing and leading
    if( sscanf(line,"%i",&fEdgeDet) != 1 )
      { PrintError(line,"<Parse edge detection mode read>"); }
    break;
  case EV1190_Window:
    // Time window, ofset and width
    if( sscanf(line,"%i%i",&fWindowOffset,&fWindow) != 2 )
      { PrintError(line,"<Parse window offset and width read>"); }
    break;
  default:
    // default try commands of TDAQmodule
    TVMEmodule::SetConfig(line, key);
    break;
  }
}

//------------------------------------------------------------------------------
void TVME_V1190::ReadIRQ( void** outBuffer )
{
  // Read and decode the memory of the V1190
  // Store the decoded ADC index and value in *outBuffer
  // Errors in contents of header word considered fatal
  // This version decodes a single event and then resets the data buffer
  //
  UShort_t status = Read(EV1190_Status);   // read status reg. for data avail
  if( !(status & 0x1) ) return;            // nothing there, return
  Int_t ich = EV1190_OutBuff;              // start register (data buffer)
  Int_t imax = ich + 1024;                 // max posible data extent
  UInt_t datum = Read(ich++);              // data header
  if( !(datum & 0x40000000) ){   // check its a header word
    ErrorStore(outBuffer, EErrDataFormat);
    ResetData();
    return;
  }
  //
  Int_t adcVal, adcIndex;                // adc value and index
  for(Int_t i=ich; i<imax; i++ ){
    datum = Read(i);
    if( (datum & 0x08000000) ) continue;    // TDC header
    if( (datum & 0x10000000) ) continue;    // TDC trailer
    if( (datum & 0x20000000) ) continue;    // TDC error
    if( (datum & 0x80000000) ) break;       // global trailer
    adcVal = datum & 0x7ffff;               // ADC value
    adcIndex = (datum & 0x3f80000) >> 19;   // ADC subaddress
    adcIndex += fBaseIndex;                 // index offset
    if( fEdgeDet == 0 ){               // pair mode
      UShort_t val1 = adcVal & 0xfff;
      UShort_t val2 = adcVal>>12;
      ADCStore(outBuffer,val1,adcIndex);
      ADCStore(outBuffer,val2,adcIndex);
    }
    else ADCStore( outBuffer, adcVal, adcIndex );// store values
  } 
  if( (datum & 0x7000000) ){   // Check trailer word
    ErrorStore(outBuffer, EErrDataFormat);
  }
  ResetData();
  return;
}

//-------------------------------------------------------------------------
void TVME_V1190::PostInit( )
{
  // Check if any general register initialisation has been performed
  // Setup data memory, threshold registers and status registers.
  // Check if ID read from internal ROM matches the expected value
  // If OK carry out default init to perform write initialisation of registers
  if( fIsInit ) return;            // Init already done?
  InitReg( V1190reg );
  TVMEmodule::PostInit();
  // Specific stuff
  // Hit detect mode...pair leading/trailing edge
  WrtMicro(EM1190_SetDetection, 0, &fEdgeDet);
  WrtMicro(EM1190_SetWinWidth, 0, &fWindow);
  WrtMicro(EM1190_SetWinOff, 0, &fWindowOffset);
  // Acquisition mode
  if( fIsContStore ) WrtMicro(EM1190_ContStore);
  else WrtMicro(EM1190_TrigMatch);
  WrtMicro(EM1190_EnSubTrig);
  return;
}

//-------------------------------------------------------------------------
Bool_t TVME_V1190::CheckHardID( )
{
  // Read module ID from internal ROM
  Int_t id = Read(EV1190_ID0) & 0xff;
  id |= (Read(EV1190_ID1) & 0xff) << 8;
  id |= (Read(EV1190_ID2) & 0xff) << 16;
  fprintf(fLogStream,"CAEN V1190 hardware ID Read: %x  Expected: %x\n",
	  id,fHardID);
  if( id == fHardID ) return kTRUE;
  else
    PrintError("","<CAEN V1190 hardware ID read error>",EErrFatal);
  return kFALSE;
}

//----------------------------------------------------------------------------
void TVME_V1190::ReadMicro( Int_t command, Int_t object, Int_t* data )
{ 
  // Read from the internal micro controller
  UInt_t opcode = (command&0xffff) | (object&0xff);
  UInt_t muHShake;
  muHShake = Read(EV1190_MicroHShake);
  if( !(muHShake & EV1190_MROK) )
    PrintError("","<Micro controller not ready for read>",EErrFatal);
  Write(EV1190_Micro, opcode);
  Int_t* dat = data;
  for(;;){
    *dat = Read(EV1190_Micro);
    dat++;
    muHShake = Read(EV1190_MicroHShake);
    if( (muHShake & EV1190_MROK) ) return;
  }
}

//----------------------------------------------------------------------------
void TVME_V1190::WrtMicro( Int_t command, Int_t object, Int_t* data )
{ 
  // Write operation to the internal micro controller
  //
  UInt_t opcode = (command&0xffff) | (object&0xff);
  UInt_t muHShake;
  for(Int_t i=0; i<EV1190_MicroTimeout; i++){
    muHShake = Read(EV1190_MicroHShake);
    if( (muHShake & EV1190_MWOK) ) break;
    usleep(1);
    if(i == EV1190_MicroTimeout)
      PrintError("","<Micro controller not ready for write>",EErrFatal);
  }
  Write(EV1190_Micro, opcode);
  if( !data ) return;
  for(Int_t i=0; i<EV1190_MicroTimeout; i++){
    muHShake = Read(EV1190_MicroHShake);
    if( (muHShake & EV1190_MWOK) ) break;
    usleep(1);
    if(i == EV1190_MicroTimeout)
      PrintError("","<Micro controller not ready for write>",EErrFatal);
  }
  Write(EV1190_Micro,data);
}

ClassImp(TVME_V1190)
