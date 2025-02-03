docker build -t bfes_image .

docker run -it bfes_image

mvn package -f /app/build/../ApacheFlink_x_Java/pom.xml

/opt/flink/bin/start-cluster.sh

/opt/flink/bin/flink run /app/ApacheFlink_x_Java/target/FlinkBFES-job-1.0-SNAPSHOT-jar-with-dependencies.jar

