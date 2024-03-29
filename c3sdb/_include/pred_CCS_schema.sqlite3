-- pred_CCS_schema.sql
-- Dylan H. Ross
-- 2019/07/17
--
--      defines the structures of the C3S.db Sqlite3 database


-- master table in which to combine all of the datasets
CREATE TABLE predicted (
    -- global unique integer identifier (same as in master, increment past highest value in master for new predictions)
    g_id INTEGER UNIQUE NOT NULL,
    -- compound name
    name TEXT NOT NULL,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- m/z and CCS (predicted)
    mz REAL NOT NULL,
    pred_ccs REAL NOT NULL,
    -- neutral smiles structure
    smi TEXT NOT NULL,
    -- class label from untargeted classification step (if used)
    class_label INTEGER,
    -- CCS error relative to reference value in master table (if available)
    pred_error REAL,
    -- timestamp (YYMMDDHHmmss) of when prediction was added to the database, NULL for original reference data
    t_stamp INTEGER
);
