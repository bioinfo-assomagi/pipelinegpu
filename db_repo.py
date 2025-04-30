import pandas as pd
from psycopg2 import sql as psycosql

from DBContext import DBContext


class DBRepo:

    def __init__(self, dbc = None):
        self.db = dbc or DBContext("dummy_path")

    def get_spawns(self, parent_id):
        sql = "SELECT riferimento as parent, sample as spawn \
                FROM acept_sample \
                WHERE sample != riferimento and riferimento = %s"
        
        results = self.db.run_query(sql, (str(parent_id), ))


        results_df = pd.DataFrame(results, columns=["parent", "spawn"])
        
        return results_df
    

class DBRepoSQLite:

    def __init__(self):
        pass