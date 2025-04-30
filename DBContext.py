import sqlite3
import psycopg2
from psycopg2 import sql as psycosql
import pandas as pd
import config
import threading

from typing import Any
from os import listdir, system


class DBContext: # TODO: if it will hold the current state of the db while the pipeline is running, a singleton may be useful
        
        def __init__(self, db_path) -> None:
            self.db_path = db_path
            
        # """ Don't allow reads, the client shouldn't bother where DBContext is pointing. """
        # def setup_database(self, db_dest):
            
        #     # TODO: add safety checks, we need to make sure we are not reading any strange db

        #     #system("rsync -avz -e ssh {} {}".format(self.server_db_mapping()[db_dest], config.DB_PATH)) #TODO: resync should be done somewhere else, where threads cannot reach (DBContext will not be used by threads right now)

        #     # if db_dest == 'b':
        #     #     self.db_path = config.DB_EUREGIO
        #     # elif db_dest == 'r':
        #     #     self.db_path = config.DB_MAGI
        #     # elif db_dest == 'z':
        #     #     self.db_path = config.DB_RICERCA

        
        # def get_disease(self, samples : list):
        #     print("DB_PATH={}".format(self.db_path))
        #     #conn = sqlite3.connect(self.db_path)

        #     conn = psycopg2.connect(dbname="limsmagidb", user="bioinfo", password="password", host="192.168.1.87", port="5432")

        #     cursor = conn.cursor()  

        #     sql = "SELECT acept_sample.sample, acept_pannelli.pannello, acept_malattia.malattia, acept_geni.gene \
        #         from acept_sample \
        #         left join acept_malattia on acept_sample.fenotipo_id = acept_malattia.id \
        #         left join acept_malattia_gene_list on acept_malattia.id = acept_malattia_gene_list.malattia_id \
        #         left join acept_geni on acept_malattia_gene_list.geni_id = acept_geni.id \
        #         left join acept_pannelli on acept_sample.panel_id = acept_pannelli.id \
        #         where acept_sample.sample in ({}) \
        #         order by sample ASC, pannello ASC, malattia ASC, gene ASC".format(','.join(['?']*len(samples))) 

        #     # results = cursor.execute(sql, samples)
        #     cursor.execute(sql, samples)
        #     results = cursor.fetchall()
        #     results_df = pd.DataFrame(results, columns=['sample', 'panel', 'malattia', 'gene'])
        #     conn.close()
        #     return results_df

        # def get_sample_familiari(self, samples):
        #     #conn = sqlite3.connect(self.db_path)
        #     conn = psycopg2.connect(dbname="limsmagidb", user="bioinfo", password="password", host="192.168.1.87", port="5432")

        #     cursor = conn.cursor()

        #     sql = "SELECT sample_id, riferimento, familiarita, familiarita_relativa, select_nucleofamiliare, AFFETTO, sesso \
        #         FROM access_accettazione \
        #         WHERE sample_id in ({})".format(','.join(['?']*len(samples)))

        #     #results = cursor.execute(sql, samples)
        #     cursor.exectue(sql, samples)
        #     results = cursor.fetchall()
        #     results_df = pd.DataFrame(results, columns=['Sample Id','Riferimento', 'Familiarita', 'Familiarita Relativa', 'Nucleo Familiare', 'Affeto', 'Sesso'])
        #     conn.close()

        #     return results_df

        def get_disease(self, samples : list):
            #print("DB_PATH={}".format(self.db_path))
            #conn = sqlite3.connect(self.db_path)

            conn = psycopg2.connect(dbname="limsmagidb", user="bioinfo", password="password", host="192.168.1.87", port="5432")

            cursor = conn.cursor()  

            sql_query = psycosql.SQL("SELECT acept_sample.sample, acept_pannelli.pannello, acept_malattia.malattia, acept_geni.gene \
                from acept_sample \
                left join acept_malattia on acept_sample.fenotipo_id = acept_malattia.id \
                left join acept_malattia_gene_list on acept_malattia.id = acept_malattia_gene_list.malattia_id \
                left join acept_geni on acept_malattia_gene_list.geni_id = acept_geni.id \
                left join acept_pannelli on acept_sample.panel_id = acept_pannelli.id \
                where acept_sample.sample in ({}) \
                order by sample ASC, pannello ASC, malattia ASC, gene ASC").format(psycosql.SQL(',').join(map(psycosql.Literal, samples))) 

            #results = cursor.execute(sql, samples)
            cursor.execute(sql_query, samples)
            results = cursor.fetchall()
            results_df = pd.DataFrame(results, columns=['sample', 'panel', 'malattia', 'gene'])
            conn.close()
            return results_df
        
        def get_sample_familiari(self, samples : list):

            conn = psycopg2.connect(dbname="limsmagidb", user="bioinfo", password="password", host="192.168.1.87", port="5432")
            cursor = conn.cursor()  

            sql_query = psycosql.SQL("SELECT sample_id, riferimento, familiarita, familiarita_relativa, select_nucleofamiliare, \"AFFETTO\", sesso \
                FROM access_accettazione \
                WHERE sample_id in ({})").format(psycosql.SQL(',').join(map(psycosql.Literal, samples))) 

            cursor.execute(sql_query, samples)
            results = cursor.fetchall()
            results_df = pd.DataFrame(results, columns=['Sample Id','Riferimento', 'Familiarita', 'Familiarita Relativa', 'Nucleo Familiare', 'Affeto', 'Sesso'])
            conn.close()

            return results_df
        
        
        def run_query(self, query, args):
            conn = psycopg2.connect(dbname="limsmagidb", user="bioinfo", password="password", host="192.168.1.87", port="5432")
            cursor = conn.cursor()

            # sql_query = psycosql.SQL(query)
            cursor.execute(query, args)
            results = cursor.fetchall()

            conn.close()

            return results
        

        """ Be careful, these are in sqlite; don't use them to access postgres. """

        def get_disease_sqlite(self, samples : list):
            print("DB_PATH={}".format(self.db_path))
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()  

            sql = "SELECT acept_sample.sample, acept_pannelli.pannello, acept_malattia.malattia, acept_geni.gene \
                from acept_sample \
                left join acept_malattia on acept_sample.fenotipo_id = acept_malattia.id \
                left join acept_malattia_gene_list on acept_malattia.id = acept_malattia_gene_list.malattia_id \
                left join acept_geni on acept_malattia_gene_list.geni_id = acept_geni.id \
                left join acept_pannelli on acept_sample.panel_id = acept_pannelli.id \
                where acept_sample.sample in ({}) \
                order by sample ASC, pannello ASC, malattia ASC, gene ASC".format(','.join(['?']*len(samples))) 

            results = cursor.execute(sql, samples)
            results_df = pd.DataFrame(results, columns=['sample', 'panel', 'malattia', 'gene'])
            conn.close()
            return results_df

        
        
        def get_sample_familiari_sqlite(self, samples):
            conn = sqlite3.connect(self.db_path)
            cursor = conn.cursor()

            sql = "SELECT sample_id, riferimento, familiarita, familiarita_relativa, select_nucleofamiliare, AFFETTO, sesso \
                FROM access_accettazione \
                WHERE sample_id in ({})".format(','.join(['?']*len(samples)))

            results = cursor.execute(sql, samples)
            results_df = pd.DataFrame(results, columns=['Sample Id','Riferimento', 'Familiarita', 'Familiarita Relativa', 'Nucleo Familiare', 'Affeto', 'Sesso'])
            conn.close()

            return results_df
        
        
        

# class DBContext:

#     class __DBContext:

#         def __init__(self) -> None:
#             self.database = None

#     _instance = None
#     _lock = threading.Lock()

#     def __init__(self):
#         if not self._instance:
#             with self._lock: # TODO: only acquire lock when needed to, i.e. when there is no instance already instatiated, when there is an _instance, threads can call this at the same time
#                 if not self._instance:
#                     self._instance = self.__DBContext()
        

#     # Delegate calls, be it method or properties, to __DBContext; since everything that a DBContext should do or contain, will be inside __DBContext
#     # DBContext is just a wrapper to implement the singleton pattern
#     def __getattr__(self, __name: str) -> Any:
#         return getattr(self._instance, __name) # e.g. if i want to call the read() method by using DBContext.read(), it will call the read() located inside __DBContext