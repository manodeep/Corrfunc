#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

//Taken from http://stackoverflow.com/questions/744766/how-to-compare-ends-of-strings-in-c
int EndsWith(const char *str, const char *suffix)
{
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

int main(int argc, char **argv)
{    
    (void) argc, (void) argv;
    if(argc != 3) {
        fprintf(stderr,"Error: Usage `%s' <text file with weights> <fastfood file with positions/velocities to append to> \n", argv[0]);
        return EXIT_FAILURE;
    }
    char *weights_file = argv[1];
    char *ff_file = argv[2];
    if(EndsWith(ff_file, ".ff") != 1) {
        fprintf(stderr,"Error: argv[2] = `%s' should be a fast-food file and end with suffix '.ff'\n",ff_file);
        return EXIT_FAILURE;
    }
    
    FILE *out = fopen(ff_file, "r+");
    if(out == NULL) {
        perror(NULL);
        return EXIT_FAILURE;
    }
    
    long endpos=0;
    if(fseek(out, 0L, SEEK_END) >= 0) {
        endpos = ftello(out);
        rewind(out);
    } else {
        fclose(out);
        perror(NULL);
        return EXIT_FAILURE;
    }

    int idat[5];
    int dummy;
    size_t nread = fread(&dummy, sizeof(dummy), 1, out);
    if(nread != 1 || dummy != 20) {
        fclose(out);
        perror(NULL);
        fprintf(stderr,"Error: Incorrect format for fast-food fortran binary. Expected to read bytes showing data for a header of 5 integers\n"
                "Read = %zu items with content = %d\n", nread, dummy);
        return EXIT_FAILURE;
    }
    nread = fread(&(idat[0]), sizeof(idat[0]), 5, out);
    if(nread != 5 || idat[1] <= 0) {
        fclose(out);
        perror(NULL);
        fprintf(stderr,"Error: Expected to read header with 5 integers. Read = %zu items with second element = %d as the number of galaxies contained (should be positive).\n",
                nread, idat[1]);
        return EXIT_FAILURE;
    }
    
    int Ngal = idat[1];

    /* We could seek to EOF and then write the weights there. However,
       re-running the code would keep adding the weights column. Therefore, 
       better to fseek to the end of the positions array and then start writing there. 
       
       However, the positions themselves could be written in double or float. 
       So, we need to first figure out how much to skip by reading in the 
       padding bytes for the positions.
     */
    long bytes_to_seek = 4 + /* end padding bytes for idat[5]*/
        (4 + 36 + 4) + /* front bytes, fdat[9], end padding bytes */
        (4 + 4 + 4);  /* front bytes, znow, end padding bytes */
        

    /* Seek to the front padding bytes for positions */
    if(fseek(out, bytes_to_seek, SEEK_CUR) < 0) {
        fclose(out);
        perror(NULL);
        return EXIT_FAILURE;
    }

    /* Check the precision of the position arrays */
    /* Read the front padding bytes */
    nread = fread(&dummy, sizeof(dummy), 1, out);
    if( ! (nread == 1 && (dummy == 4*Ngal || dummy == 8*Ngal))) {
        fclose(out);
        perror(NULL);
        fprintf(stderr,"Error: Expected to read one item of padding bytes, read = %zu items. Positions can be either float (4 bytes) or doubles (8 bytes) - found %d instead.\n"
                "Ngal = %d padding bytes = %d\n",
                nread, dummy/Ngal, Ngal, dummy);
        return EXIT_FAILURE;
    }
    long curr_pos = ftell(out);
    const int precision = dummy/Ngal;

    
    bytes_to_seek = dummy + 4 + /* we know how many bytes already, end padding bytes*/
        + (4 + dummy + 4) + /* front, 2nd axis positions, end padding bytes */
        + (4 + dummy + 4); /* front, 3rd axis positions, end padding bytes */

    if(truncate(ff_file, curr_pos + bytes_to_seek) < 0) {
        perror(NULL);
        fprintf(stderr,"Error: Could not truncate file to correct size = %ld \n", curr_pos + bytes_to_seek);
        fclose(out);
        return EXIT_FAILURE;
    }
    
#if 0    
    /* gals_Mr19.ff actually contained velocities (I think, could have been 3 fields of halomass, central/sat, haloid) */
    const long bytes_for_vel_fields = (4 + dummy + 4) + /* front, 1st axis velocities, end padding bytes */
        + (4 + dummy + 4) + /* front, 2nd axis velocities, end padding bytes */
        + (4 + dummy + 4); /* front, 3rd axis velocities, end padding bytes */
    
    const int skip_velocities = (curr_pos + bytes_to_seek + bytes_for_vel_fields) <= endpos;/* if weights has been written before, then the "<" will hold. */
    if(skip_velocities == 1) {
        fprintf(stderr,"Skipping velocity fields\n");
        bytes_to_seek += bytes_for_vel_fields;
    }
#endif    

    /* Seek to the end of positions and velocities */
    if(fseek(out, bytes_to_seek, SEEK_CUR) < 0) {
        fclose(out);
        perror(NULL);
        return EXIT_FAILURE;
    }
    long pos = ftell(out);
    fprintf(stderr,"Current file position = %ld end pos = %ld\n", pos, endpos);
    

    FILE *in = fopen(weights_file, "r");
    if(in == NULL){
        fclose(out);
        perror(NULL);
        return EXIT_FAILURE;
    }

    char buf[128];
    float flt_weight;
    double dbl_weight;
    int nweights=0;
    dummy=Ngal*precision;
    size_t nwritten = fwrite(&dummy, sizeof(dummy), 1, out);
    if(nwritten != 1) {
        fclose(out);
        perror(NULL);
        return EXIT_FAILURE;
    }
    while(fgets(buf, 128, in) != NULL) {
        /* fprintf(stderr,"buf = %s\n",buf); */
        dbl_weight = atof(buf);
        flt_weight = (float) dbl_weight;
        if(precision == 4) {
            nwritten = fwrite(&flt_weight, sizeof(flt_weight), 1, out);
        } else {
            nwritten = fwrite(&dbl_weight, sizeof(dbl_weight), 1, out);
        }
        if(nwritten != 1) {
            fclose(out);
            perror(NULL);
        }
        nweights++;
    }
    nwritten = fwrite(&dummy, sizeof(dummy), 1, out);
    if(nwritten != 1) {
        fclose(out);
        perror(NULL);
    }
    
    fclose(in);
    fclose(out);

    fprintf(stderr,"Success: Wrote %d weights into file. Extra size = %d \n",nweights, dummy);
    return EXIT_SUCCESS;
}
