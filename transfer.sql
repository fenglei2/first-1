USE prognosticsys;

load data local infile 
"C:/Users/nian/Desktop/BreastDataSet.txt"
into table breastdataset
(Grade, size, nodes, ER, agedx, race, vs, time, delta);   
DELETE FROM breastdataset limit 1;
alter table breastdataset drop column vs;  

UPDATE breastdataset SET delta = 0 WHERE delta!=46;
UPDATE breastdataset SET delta = 1 WHERE delta=46;

UPDATE breastdataset SET agedx = 1 WHERE agedx<=50;
UPDATE breastdataset SET agedx = 2 WHERE agedx>50;

DELETE FROM breastdataset WHERE ER!=1 AND ER!=2;

DELETE FROM breastdataset WHERE nodes>20;
UPDATE breastdataset SET nodes = 2 WHERE nodes>=1 AND nodes<=3;
UPDATE breastdataset SET nodes = 3 WHERE nodes>=4 AND nodes<=10;
UPDATE breastdataset SET nodes = 4 WHERE nodes>10 AND nodes<=20;
UPDATE breastdataset SET nodes = 1 WHERE nodes=0;


DELETE FROM breastdataset WHERE size>100;
UPDATE breastdataset SET size = 1 WHERE size<=20;
UPDATE breastdataset SET size = 2 WHERE size>20 AND size<=50;
UPDATE breastdataset SET size = 3 WHERE size>50 AND size<=100;

DELETE FROM breastdataset WHERE Grade!=1 AND Grade!=2 AND Grade!=3 AND Grade!=4;
UPDATE breastdataset SET Grade = 3 WHERE Grade=4;




