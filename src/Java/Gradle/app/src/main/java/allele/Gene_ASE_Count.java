package allele;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;


public class Gene_ASE_Count
{
    protected int Good; //Keeps track of the number of reads that support a single allele (want to be large)
    protected int Bad; //Keeps track of the number of reads that support multiple alleles (want to be small)
    protected FileWriter savFil; //The file the output will be saved to
    protected File bamFile; //The bamFile to be processed
    protected String gene_tag; //The tag in the bam file used to ID the gene or other region of interest (peak, etc)
    protected HashMap<String, String> UMIMap; //Dictionary mapping each cbc/UMI/gene triplet to the allele it corresponds to
    protected HashMap<String,Integer> counts; //Dictionary containing information on the number of UMI supporting each cbc/gene/allele triplet


    //Generates the initial object and parses arguments
    public Gene_ASE_Count(String args[]) throws Exception
    {
        //ArrayList<String[]> ret=new ArrayList<String[]>();
        this.bamFile = new File(args[0]); 
        //File vcfFile=new File(args[1]);
        //int ind=Integer.parseInt(args[2]);
        String saveFile=args[1];

        this.gene_tag="GN";
        if(args.length==3)
        {
            this.gene_tag=args[2];
        }

        this.savFil = new FileWriter(saveFile); 
    }


    //Reads through all reads one by one and processes them
    public void IterateOverAll() throws Exception
    {

        print("Read in long read data!");
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(this.bamFile);
        SAMRecordIterator r = sr.iterator();
        int iter=0;
        this.Good=0;
        this.Bad=0;


        this.UMIMap=new HashMap<String, String>(); 

        while(r.hasNext()) {
            iter=iter+1;

            if(iter%10000==0)
            {
                System.out.println(this.Good);
                System.out.println(this.Bad);
                System.out.println((float)this.Good/(float)(this.Good+this.Bad));
                System.out.println(iter);
            }

            SAMRecord read=r.next();//current read

            ProcessRead(read);

        }

        System.out.println("Get Allele Counts");

        UMIMAPtoCounts();

        System.out.println("Save!");

        saveCounts();
       
        System.out.println("Clean up");

        r.close();
        sr.close();        
        this.savFil.close();

    }

    //Saves the UMI counts 
    public void saveCounts() throws Exception
    {
        Iterator<String> it_cnt=this.counts.keySet().iterator();

        while (it_cnt.hasNext()) {
            String key=it_cnt.next();
            Integer val=this.counts.get(key);

            String res=key+" "+Integer.toString(val);
            this.savFil.write(res+ System.lineSeparator());

        }

    }

    //Takes the internal UMIMAP object and gets UMI count information from it
    public void UMIMAPtoCounts()
    {
        this.counts=new HashMap<String,Integer>();
        Iterator<String> it=this.UMIMap.keySet().iterator();


        while (it.hasNext()) {
            String key=it.next();
            String val=this.UMIMap.get(key);
            String[] split=key.split(" ");
            String cbc=split[1];
            String gene=split[2];
            //String allele=split[3];
            String res=cbc+" "+gene+" "+val;
            //savFil.write(res+ System.lineSeparator());
            if(!this.counts.containsKey(res))
            {
                this.counts.put(res,0);
            }

            this.counts.put(res,this.counts.get(res)+1);

        }
    }


    //processes and individual read to get allele information and to filter out certain reads
    public void ProcessRead(SAMRecord read)
    {

        //Get needed info from read
        String cbc=read.getStringAttribute("CB");
        String umi=read.getStringAttribute("UB");
        String gene=read.getStringAttribute(gene_tag);

        int numb=0;//number alignments for read

        try{
            numb=read.getIntegerAttribute("NH");
        }
        catch(Exception e)
        {
            return;
        }

        int WASP=0; //Score according to WASP
        try{
            WASP=read.getIntegerAttribute("vW");
        }
        catch(Exception e)
        {
            return;
        }


        //Check passes filters, if not skip read
        if(WASP!=1)
        {
            return;
        }

        if(gene==null || cbc==null || umi==null)
        {
            return;
        }

        if(numb>1)
        {
            return;
        }

        //Get allele specific information
        int All1=0;//Number SNPs supporting allele 1
        int All2=0; //Number SNPs supporting allele 2
        int Tot=0; //Total number of SNPs

        byte[] Alls=read.getSignedByteArrayAttribute("vA"); //read in

        Tot=Alls.length;
        if(Tot<1)
        {
            return;
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


        if(Tot<1){return;}
        
        float rat=(float)Math.max(All1,All2)/(float)Tot;//percent of reads supporting allele

        if(rat>.95)
        {
            this.Good=this.Good+1;
        }
        else
        {
            this.Bad=this.Bad+1;
            return;
        }

        //Adds allele information to UMIMap
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

        if(this.UMIMap.containsKey(res))
        {
            String cur=this.UMIMap.get(res);
            if(!cur.equals(val))
            {
                this.UMIMap.put(res,"Ambig");
            }
        }
        else
        {
            this.UMIMap.put(res,val);
        }


    }

    //Used to print things in a less verbose manner
    public static void print(String toPrint)
    {

        System.out.println(toPrint);

    }

}

