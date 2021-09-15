package Phasing_Illumina;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;


public class PhaseCount_WASP
{
public static void main(String args[]) throws Exception
{
ArrayList<String[]> ret=new ArrayList<String[]>();
File bamFile = new File(args[0]); 
File vcfFile=new File(args[1]);
int ind=Integer.parseInt(args[2]);
String saveFile=args[3];

FileWriter savFil = new FileWriter(saveFile); 


print("Read in long read data!");
SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
SAMRecordIterator r = sr.iterator();
int iter=0;
int Good=0;
int Bad=0;


HashMap<String, String> UMIMap=new HashMap<String, String>(); 

while(r.hasNext()) {
iter=iter+1;

if(iter%10000==0)
{
System.out.println(Good);
System.out.println(Bad);
System.out.println((float)Good/(float)(Good+Bad));
System.out.println(iter);
}

SAMRecord read=r.next();

String seq=read.getReadString();

String cbc=read.getStringAttribute("CB");
String umi=read.getStringAttribute("UB");
String gene=read.getStringAttribute("GN");

int numb=0;
try{
numb=read.getIntegerAttribute("NH");
}
catch(Exception e)
{
continue;
}

int WASP=0;
try{
WASP=read.getIntegerAttribute("vW");
}
catch(Exception e)
{
continue;
}

if(WASP!=1)
{
continue;
}

if(gene==null || cbc==null || umi==null)
{
continue;
}

if(numb>1)
{
continue;
}

int All1=0;
int All2=0;
int Tot=0;

///////////////////
///New stuff added here!
/////////////////////
byte[] Alls=read.getSignedByteArrayAttribute("vA"); //read in

Tot=Alls.length;
if(Tot<1)
{
continue;
}

for(int i=0;i<Tot;i++)
{
byte val=Alls[i];
if(val==1)
{
All1=All1+1;
}
if(val==2)
{
All2=All2+1;
}
}


if(Tot<1){continue;}
//if(Tot!=All1+All2){print("Yuck!");continue;} //If any do not match throws out read
float rat=(float)Math.max(All1,All2)/(float)Tot;

if(rat>.95)
{
Good=Good+1;
}
else
{
Bad=Bad+1;
continue;
}



String res=umi+" "+cbc+" "+gene;
String val="None";
if(All1>All2)
{
val="All1";
}
if(All2>All1)
{
val="All2";
}

if(UMIMap.containsKey(res))
{
String cur=UMIMap.get(res);
if(!cur.equals(val))
{
UMIMap.put(res,"Ambig");
}
}
else
{
UMIMap.put(res,val);
}




}

System.out.println("Get Allele Counts");

HashMap<String,Integer> counts=new HashMap<String,Integer>();
Iterator<String> it=UMIMap.keySet().iterator();


while (it.hasNext()) {
String key=it.next();
String val=UMIMap.get(key);
String[] split=key.split(" ");
String cbc=split[1];
String gene=split[2];
//String allele=split[3];
String res=cbc+" "+gene+" "+val;
//savFil.write(res+ System.lineSeparator());
if(!counts.containsKey(res))
{
counts.put(res,0);
}

counts.put(res,counts.get(res)+1);

}

System.out.println("Save!");


Iterator<String> it_cnt=counts.keySet().iterator();

while (it_cnt.hasNext()) {
String key=it_cnt.next();
Integer val=counts.get(key);

String res=key+" "+Integer.toString(val);
savFil.write(res+ System.lineSeparator());

}

r.close();
sr.close();
savFil.close();

}




public static void print(String toPrint)
{

System.out.println(toPrint);

}

}

