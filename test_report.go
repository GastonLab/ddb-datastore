/*
This is a testing program in GO designed to generate clinical genomics reports.
This is expected to be more performant than the python scripts and to make better
use of parallele processing in filtering and categorizing variants from the database.
*/

package main

import (
  "fmt"
  "log"

  "github.com/gocql/gocql"
)

func main () {
  // Connect to the cassandra cluster
  cov_cluster := gocql.NewCluster("142.239.155.181", "142.239.155.182",
    "142.239.155.183", "142.239.155.184")
  cov_cluster.Keyspace = "coveragestore"
  cov_cluster.Consistency = gocql.Quorum
	cov_session, _ := cov_cluster.CreateSession()
	defer cov_session.Close()

  var_cluster := gocql.NewCluster("142.239.155.181", "142.239.155.182",
    "142.239.155.183", "142.239.155.184")
  var_cluster.Keyspace = "variantstore"
  var_cluster.Consistency = gocql.Quorum
	var_session, _ := var_cluster.CreateSession()
	defer var_session.Close()

  // list all tweets (Copied from Tutorial, not adapted yet)
	iter := session.Query(`SELECT id, text FROM tweet WHERE timeline = ?`, "me").Iter()
	for iter.Scan(&id, &text) {
		fmt.Println("Tweet:", id, text)
	}
	if err := iter.Close(); err != nil {
		log.Fatal(err)
	}
}
