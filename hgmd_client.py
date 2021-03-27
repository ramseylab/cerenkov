"""
Created on Aug 3, 2016

@author: ramseylab
"""

import pandas
from sqlalchemy import create_engine
from sqlalchemy.engine.url import URL
from sqlalchemy.pool import NullPool


class HgmdClient:
    # Make sure you have created such a user in MySQL
    # For MySQL 5.6
    #   GRANT SELECT PRIVILEGES ON hg19.* To 'bud'@'localhost' IDENTIFIED BY 'earth';
    #   GRANT SELECT PRIVILEGES ON hgmd_pro.* To 'bud'@'localhost' IDENTIFIED BY 'earth';
    # For MySQL 5.7
    #   CREATE USER 'bud'@'localhost' IDENTIFIED BY 'earth';
    #   GRANT SELECT ON hg19.* TO 'bud'@'localhost';
    #   GRANT SELECT ON hgmd_pro.* TO 'bud'@'localhost';
    #   FLUSH PRIVILEGES;
    __db_url = dict(
        drivername='mysql+pymysql',
        host='localhost',
        port='3306',
        username='bud',
        password='earth',
        database='hgmd_pro',
        query={'charset': 'utf8'}
    )

    def __init__(self):
        # db = create_engine('mysql://bud:earth@localhost:3306/hg19') # require module `MySQLdb`
        #   default dialect is 'mysql+mysql-python'
        #   `MySQLdb` is a fork of MySQL-python with added support for Python 3
        #   See http://docs.sqlalchemy.org/en/latest/core/engines.html#mysql

        # db = create_engine('mysql+pymysql://bud:earth@localhost:3306/hg19') # require module `PyMySQL`

        # For `poolclass`, see http://stackoverflow.com/a/8705750

        self.db = create_engine(URL(**HgmdClient.__db_url), poolclass=NullPool)
        self.conn = self.db.connect()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()
        self.db.dispose()

    def select_rsid(self):
        # Use '%%' to escape '%' in python
        query = '''
                SELECT DISTINCT dbSNP as name
                FROM allmut
                WHERE base IN ('R')
                      AND disease NOT LIKE '%%cancer%%'
                AND dbSNP IS NOT NULL
                '''

        rows = self.conn.execute(query)

        return [x[0] for x in rows.fetchall()]
