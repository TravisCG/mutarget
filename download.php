<?php
  $suffix = $_GET["name"];
  if (substr($suffix, 0, 5) == "sites"){
    $filename = substr($suffix, 27);
  }
  else {
    $filename = "/tmp/MUT" . $suffix;
  }
  header('Content-Type: text/plain');
  $f = fopen($filename . ".tsv", "r");
  while($line = fgets($f)){
	  echo $line;
  }
  fclose($f);
?>
