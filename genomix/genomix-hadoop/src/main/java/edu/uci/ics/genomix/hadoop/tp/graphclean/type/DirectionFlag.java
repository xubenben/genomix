package edu.uci.ics.genomix.hadoop.tp.graphclean.type;

import edu.uci.ics.genomix.hadoop.tp.graphclean.type.WholeStateInformation.State;

public class DirectionFlag extends State {
    public static final byte DIR_NO = 0b000 << 0;
    public static final byte DIR_FF = 0b001 << 0;
    public static final byte DIR_FR = 0b010 << 0;
    public static final byte DIR_RF = 0b011 << 0;
    public static final byte DIR_RR = 0b100 << 0;
    public static final byte DIR_MASK = 0b111 << 0;
    public static final byte DIR_CLEAR = 0b1111000 << 0;
    
//    public static final byte DIR_FROM_DEADVERTEX = 0b101 << 0;
}