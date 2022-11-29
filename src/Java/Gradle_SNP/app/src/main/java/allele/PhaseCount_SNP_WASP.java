package allele;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;


public class PhaseCount_SNP_WASP
{
    public static void main(String args[]) throws Exception
    {
        SNP_ASE_Count ASECounter=new SNP_ASE_Count(args);
        ASECounter.IterateOverAll();

    }




    public static void print(String toPrint)
    {

        System.out.println(toPrint);

    }

}

