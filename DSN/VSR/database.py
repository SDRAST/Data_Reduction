"""
Functions to put VSR data in the database

The data to be written to the database are obtained with parse_VSR_log()
which returns a nest set of dictionaries that looks like this::
  {"1A":{"1W1":{data},"1W2":{data}}, "1B":{"1W1":{data},"2W1":{data}}}

What will be written to the MySql table for each subchannel is::
  "1A","1W1",...data...
  "1A","2W1",...data...
  "2A","1W1",...data...
  "2A","2W1",...data...

The MySql queries to add the necessary tables to the database are::

  -- venus.vsr_pars
  create table venus.vsr_pars (
    ID                        int unsigned auto_increment not null,
    vsr                       varchar(8) not null,
    subchannel                varchar(4) not null,
    LO                        int not null,
    bpsamp                    tinyint not null,
    bw                        int not null,
    datafile                  varchar(128) null,
    dss                       tinyint,
    frov                      decimal(20,6) not null,
    restfreq                  decimal(8,1) not null,
    sfro                      decimal(12,1),
    constraint pk_Venus_vsr_pars primary key (ID)
  ) type = InnoDB ;
  
  -- venus.vsr_times
  create table venus.vsr_times (
    ID                        int unsigned auto_increment not null,
    start                     decimal(13,1) not null,
    stop                      decimal(13,1) not null,
    vsr_pars_id               int unsigned not null,
    constraint pk_Venus_vsr_times primary key (ID)
  ) type = InnoDB ;
  
  alter table venus.vsr_times add constraint fk_vsr_pars
    foreign key (vsr_pars_id)
    references venus.vsr_pars (ID) on update cascade ;
"""
 
import Mysql as sql

def store_VSR_data(db,database,new_table,data):
  """
  Store VSR data in DSN database

  @param db : a database object

  @param new_table : boolean -
    If True, the current table is dropped and re-created.
    
  @param data : dictionary
    Data as returned from parse_VSR_log()

  @return: None
  """
  if new_table:
    response = sql.ask_db(db,"drop table if exists "+database+".vsr_times ;")
    response = sql.ask_db(db,"drop table if exists "+database+".vsr_pars ;")
    response = sql.ask_db(db,
    "create table " + database + ".vsr_pars (ID int unsigned auto_increment not null, vsr varchar(8) not null, subchannel varchar(4) not null, LO int not null, bpsamp tinyint not null, bw int not null, datafile varchar(128) null, dss tinyint, frov decimal(20,6) not null, restfreq decimal(8,1) not null, sfro decimal(12,1), constraint pk_" + database.capitalize()+"_vsr_pars primary key (ID)) ENGINE = InnoDB ;")
    response = sql.ask_db(db,"create table " + database + ".vsr_times (ID int unsigned auto_increment not null, start decimal(13,1) not null, stop decimal(13,1) not null, vsr_pars_id int unsigned not null, constraint pk_" + database.capitalize() + "_vsr_times primary key (ID)) ENGINE = InnoDB ;")
    response = sql.ask_db(db,"alter table "+database+".vsr_times add constraint fk_vsr_pars foreign key (vsr_pars_id) references "+database+".vsr_pars (ID) on update cascade ;")
  # 
  par_fields = {}
  time_fields = {}
  for vsr in data.keys():
    # Do all channels
    par_fields['vsr'] = vsr
    for subchan in data[vsr].keys():
      # Do all subchannels for each channel
      par_fields['subchannel'] = subchan
      for name in data[vsr][subchan].keys():
        # copy each parameter for each subchannel into a record
        if name != "start" and name != "stop":
          par_fields[name] = data[vsr][subchan][name]
      pars_ID = sql.insert_record(db,database+".vsr_pars", par_fields)
      print "New VSR pars record",pars_ID[0]
      # There may be multiple start/time times
      for i in range(len(data[vsr][subchan]["start"])):
        time_fields["start"] = data[vsr][subchan]["start"][i]
        time_fields["stop"] = data[vsr][subchan]["stop"][i]
        time_fields["vsr_pars_id"] = pars_ID[0]
        time_ID = sql.insert_record(db,database+".vsr_times", time_fields)
        print "New VSR times record",time_ID[0]

