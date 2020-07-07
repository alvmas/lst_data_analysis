#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 09:33:23 2019

@author: alvarom
"""

'''
Based on c code from M.Punch
'''

'''
#Define the different parameters that will be necessary in the communication

UCTS_ADDRESS= "10.10.4.3"   #GAE LAB

UCTS_SOCKET_ADDRESS= "10.10.7.250" #GAE LAB

#define CAMERA_SERVER_ADDRESS "147.96.10.41" // GAE LAB

#number of events in one UCTS packet
NB_EVENT_IN_PACKET_UCTS=20
'''


import socket
import sys
import os 
from bitarray import bitarray #IMPORTANT TO INSTALL THE MODULE (pip install bitarray)
import struct
import time
from timing_tests import execute_tests
from timing_tests import read_file


def receive_buffers(rcvBufferDescriptors,rcvBuffer):
    
    #UCTS_ADDRESS= "10.10.4.5" 
    UCTS_ADDRESS= "10.10.4.5"
    #Adress where buffers are sent from the UCTS
    #UCTS_SOCKET_ADDRESS= "10.10.7.250"
    #UCTS_SOCKET_PORT=55000 
    
    UCTS_SOCKET_ADDRESS= "10.10.7.250"
    UCTS_SOCKET_PORT=55000 
    #Adress of the socket
    SocketAddressPort=(UCTS_SOCKET_ADDRESS, UCTS_SOCKET_PORT) 
    
    # CREATE a socket using the function socket.socket
    SocketUCTS = socket.socket(family=socket.AF_INET, type=socket.SOCK_DGRAM)
    print('Socket created')
    
    #Bind the socket to the port
    SocketUCTS.bind(SocketAddressPort)
    print('Socket connected to', SocketAddressPort)
    
    #Set a timeout
    timeout=60
    SocketUCTS.settimeout(timeout)
    
    head_receive=0
    size=0
    packets_received=0
    packets_lost=0
    tail_receive=RCV_BUFFER_NUMBER-1 #End of the circular buffer
    state=True
    
    while size>=0 and state: #state=False or True
        
        #Buffer in which we want to store the information
        buffer_to_write=head_receive%RCV_BUFFER_NUMBER
        
        #Read the information sent by the UCTS
        rcvBufferDescriptors[buffer_to_write] =  SocketUCTS.recvfrom_into(rcvBuffer[buffer_to_write])
        size = rcvBufferDescriptors[buffer_to_write][0]
            

        print('Buffer received with size=',size,'bytes')
        
        if size>0 and str(rcvBufferDescriptors[buffer_to_write][1][0]) == UCTS_ADDRESS:
            if True:
                #Divide into events and extract the different bits	
                div_buffer(rcvBuffer[buffer_to_write],size)
                print('Packet received')
                packets_received+=1
                head_receive+=1
            else:
                packets_lost+=1    
                print('Packet lost')
        
      
    #Close the socket
    SocketUCTS.close()
    
   
def extract_bits(buff, offset, length): 
    #Calculate the starting byte,the ending byte 
    #and the right shift depending on the length we want to extract
    
    byte_offset_start = int(offset/8)
    byte_offset_end = int((offset+length-1)/8)
    right_shift= (byte_offset_end+1)* 8 - (offset + length)
    
    #We concatenate all the bytes in byte_bits from the byte_start to the byte_end 
    byte_bits=0
    for i in range(byte_offset_start,byte_offset_end+1):
        byte_bits = byte_bits|buff[i]
        if i < byte_offset_end:
            byte_bits = byte_bits << 8
            
    #Use a shift and a mask to obtain only the corresponding bits    
    byte_bits = (byte_bits >> right_shift) & (int('1'*length,2))

    return byte_bits


def extract_bits(buff,offset, length): 
    
    byte_offset_start = int(offset/8)
    byte_offset_end = int((offset+length-1)/8)
    left_shift = offset - byte_offset_start*8
    right_shift= (byte_offset_end+1)*8-(offset+length)
    #print(bin(buff[byte_offset_start]))
    
    if byte_offset_end == byte_offset_start:
        bits_extracted=(buff[byte_offset_start]>>right_shift) & (int('1'*length,2))
        return bits_extracted   
    else:
        byte_bits=buff[byte_offset_start]<<left_shift
        byte_bits=byte_bits<<(8-left_shift)
        
        for i in range(byte_offset_start+1, byte_offset_end):
            byte_bits=byte_bits|buff[i]
            byte_bits=byte_bits<< 8
            
            
        #Now we go for the last byte to extract from
        last_byte_bits=buff[byte_offset_end] >> right_shift
        last_byte_bits=last_byte_bits<< right_shift
        bits_extracted=byte_bits|last_byte_bits
        bits_extracted=bits_extracted >> right_shift
        
        return bits_extracted    

def writetxt(event_info,name):
     wr_file = open (name, "a")
     for i in range(0,len(event_info['eventcounter_readout'])):
	     wr_file.write(str(event_info['bch'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['eventcounter_readout'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['eventcounter_busy'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['timestamp'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['ppscounter'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['triggerType'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['telescope_pattern'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['spi'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['clk'][i])+'\n')
	     wr_file.write('\t')
	     wr_file.write(str(event_info['flg_tm_time_valid'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['flg_rst_cnt_ack'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['flg_Busy'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['flg_SPLL'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['version'][i]))
	     wr_file.write('\t')
	     wr_file.write(str(event_info['system_time'][i])+'\n')

     wr_file.close()





     
def createtxt(name):
    wr_file = open (name, "w")
    wr_file.write(time.strftime("%c")+'\n')
    wr_file.write('Bch')
    wr_file.write('\t')
    wr_file.write('Evt_readout')  
    wr_file.write('\t')
    wr_file.write('Evt_busy')  
    wr_file.write('\t')
    wr_file.write('TimeStamp')
    wr_file.write('\t')
    wr_file.write('pps') 
    wr_file.write('\t')
    wr_file.write('triggerType') 
    wr_file.write('\t')
    wr_file.write('telescope_pattern') 
    wr_file.write('\t')
    wr_file.write('spi') 
    wr_file.write('\t')
    wr_file.write('clk')
    wr_file.write('\t')
    wr_file.write('flg_tm_time_valid') 
    wr_file.write('\t')
    wr_file.write('flg_rst_cnt_ack') 
    wr_file.write('\t')
    wr_file.write('flg_Busy') 
    wr_file.write('\t')
    wr_file.write('flg_SPLL') 
    wr_file.write('\t')
    wr_file.write('Version') 
    wr_file.write('\t')
    wr_file.write('system_time' +'\n')
    wr_file.close()


    
def div_buffer(buffer,size):
    
    #Primero dividimos los datos en tail y en eventos
    events_length=size-TAIL_SIZE
    
    #Only if we have an integer number of events we read the buffer
    if events_length==0:
        print('No events in this bunch')
    else:
        if events_length%EVT_SIZE==0:
            number_events=int(events_length/EVT_SIZE)
            print(number_events)
            #For each event, divide the data into the different events + tailer
            data=(struct.unpack_from(str(size)+'B', buffer))
	   
		
            #Extract the events from the buffer and store them in packet_events
            packet_events=[]
            for i in range(0,number_events):
                events_data=[]
                for j in range(i*EVT_SIZE,(i+1)*EVT_SIZE):
            		    events_data.append(data[j])
                packet_events.append(events_data)
		
		
            #Extract the events from the buffer and store them in tailer
            tailer_data=[]
            for j in range(number_events*EVT_SIZE,len(data)):
                tailer_data.append(data[j])
            
            
            
            #Extract the tailer bits   
            bch_counter     = extract_bits(tailer_data, 0,  32)
            evt_MSB_readout = extract_bits(tailer_data, 32, 32)
            evt_MSB_busy = extract_bits(tailer_data, 64, 32)
            pps_MSB = extract_bits(tailer_data, 96, 16)
            sec_MSB = extract_bits(tailer_data, 112, 32)
            flg_tm_time_valid  = extract_bits(tailer_data, 144, 1)
            flg_rst_cnt_ack  = extract_bits(tailer_data, 145, 1)
            other_info  = extract_bits(tailer_data, 146, 6)
            version = extract_bits(tailer_data, 152, 8)
	    
            
            
            #Extract the information from the last event in the bunch
            evt_readout_final=extract_bits(packet_events[-1], 16, 8)
            evt_busy_final=extract_bits(packet_events[-1], 24, 8)
            pps_final=extract_bits(packet_events[-1], 32, 2)
            sec_final=extract_bits(packet_events[-1], 34, 2)
            

            #Extract from each packet the different bits
            for i in range(0,len(packet_events)):

                spi     = extract_bits(packet_events[i],  0, 16)
                evt_readout_LSB = extract_bits(packet_events[i], 16, 8)
                evt_busy_LSB = extract_bits(packet_events[i], 24, 8)
                pps_LSB = extract_bits(packet_events[i], 32,  2)
                sec_LSB = extract_bits(packet_events[i], 34, 2)
                flg_Busy=extract_bits(packet_events[i], 36, 1)
                flg_SPLL=extract_bits(packet_events[i], 37, 1)
                clk  = extract_bits(packet_events[i], 38, 26)
                ens     = extract_bits(packet_events[i], 64, 28)
                other_info_event     = extract_bits(packet_events[i], 92, 1)
                ons     = extract_bits(packet_events[i], 93, 3)
                

                event_readout=(((evt_MSB_readout)>>8)-(evt_readout_final<evt_readout_LSB))<<8|evt_readout_LSB
                event_busy=(((evt_MSB_busy)>>8)-(evt_busy_final<evt_busy_LSB))<<8|evt_busy_LSB
                pps=(((pps_MSB)>>2)-(pps_final<pps_LSB))<<2|pps_LSB
                sec=(((sec_MSB)>>2)-(sec_final<sec_LSB))<<2|sec_LSB
		 
                print("NUEVO EVENTO")

                #We store in a dictionary
                event_info['eventcounter_readout'].append(event_readout)
                event_info['eventcounter_busy'].append(event_busy)
                event_info['ppscounter'].append(pps)
                event_info['clk'].append(clk)
                event_info['timestamp'].append(sec*10**9+ens*8+ons)
                event_info['bch'].append(bch_counter)
                event_info['triggerType'].append(spi & 255)
                event_info['system_time'].append(int(time.time()*10**9))
		event_info['flg_tm_time_valid'].append(flg_tm_time_valid)
		event_info['flg_rst_cnt_ack'].append(flg_rst_cnt_ack)
		event_info['telescope_pattern'].append(spi>>8)
		event_info['flg_Busy'].append(flg_Busy)
		event_info['flg_SPLL'].append(flg_SPLL)
		event_info['spi'].append(spi)
		event_info['version'].append(version)


                
                #We write into a txt file the different events and print them
                #writetxt(event_dic)
                
                if True:
                    print('spi=',spi)
		    print('TriggerType=',str(spi & 255))
		    print('Telescope_pattern=',str(spi >>8))
                    print('evt_readout=',event_readout)
                    print('evt_b=',event_busy)
                    print('pps=',pps)
                    print('clk=',clk)
                    print('sec=',sec)
                    print('8ns=',ens)
                    print('1ns=',ons)
                    print('timestamp=',sec*10**9+ens*8+ons)
		    print('version=',version)
		      
        else:
            print('Does not correspond to a integer number of events')
            return
       

if __name__ == '__main__':
    


    
    #Define the Buffer Descriptor and the Buffer information
    RCV_BUFFER_NUMBER=1000
    RCV_BUFFER_SIZE= 1500
    rcvBufferDescriptors=[bytearray(RCV_BUFFER_SIZE)] * RCV_BUFFER_NUMBER
    rcvBuffer=[bytearray(RCV_BUFFER_SIZE)] * RCV_BUFFER_NUMBER
    rcvBuffer2=[bytearray(RCV_BUFFER_SIZE)] * RCV_BUFFER_NUMBER
    
    
    #Each Buffer has the same structure (24 events (syze=12 bytes)+tailer (20 bytes))
    TAIL_SIZE=20 
    EVT_SIZE=12

    #Define the lists to concatenate the different values for all events
    event_info={}
    event_info['eventcounter_readout']=[]
    event_info['eventcounter_busy']=[]
    event_info['ppscounter']=[]
    event_info['clk']=[]
    event_info['timestamp']=[]
    event_info['bch']=[]
    event_info['triggerType']=[]
    event_info['arbitraryInformation']=[]
    event_info['system_time']=[]
    event_info['flg_tm_time_valid']=[]
    event_info['flg_rst_cnt_ack']=[]
    event_info['telescope_pattern']=[]
    event_info['flg_Busy']=[]
    event_info['flg_SPLL']=[]
    event_info['spi']=[]
    event_info['version']=[]

    
    #Call the function to read the data from the UCTS
    try:
    	receive_buffers(rcvBufferDescriptors,rcvBuffer)
    except KeyboardInterrupt:
    	name='Eventos'+str(time.strftime("%d_%m_%y"))+'_'+str(time.strftime("%X"))+'.txt'
	createtxt(name)
	writetxt(event_info,name)
	event_dic=read_file(name)
	execute_tests(event_dic,True)
	exit()

