#include <sys/resource.h>
#include <stdio.h>
#include <mach/mach.h>
/*typedef struct {*/
    /*unsigned long size,resident,share,text,lib,data,dt;*/
/*} statm_t;*/

/*void read_off_memory_status(statm_t& result)*/
/*{*/
  /*unsigned long dummy;*/
  /*const char* statm_path = "/proc/self/statm";*/

  /*FILE *f = fopen(statm_path,"r");*/
  /*if(!f){*/
    /*perror(statm_path);*/
    /*abort();*/
  /*}*/
  /*if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",*/
    /*&result.size,&result.resident,&result.share,&result.text,&result.lib,&result.data,&result.dt))*/
  /*{*/
    /*perror(statm_path);*/
    /*abort();*/
  /*}*/
  /*fclose(f);*/
/*}*/

/*int main() {*/
void pmem(){
  struct rusage r_usage;
  getrusage(RUSAGE_SELF,&r_usage);
  printf("Memory usage: %ld bytes\n",r_usage.ru_maxrss);
  printf("Memory usage: %ld kilobytes\n",r_usage.ru_maxrss/1024);
  printf("Memory usage: %ld megabytes\n",r_usage.ru_maxrss/1024/1024);

  /*return 0;*/
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        printf("can't access");
    printf(" mem %ld bytes \n",(size_t)info.resident_size);

        /*return (size_t)0L;      [> Can't access? <]*/
    /*return (size_t)info.resident_size;*/
}


