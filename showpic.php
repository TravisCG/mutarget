<?php
  $suffix = $_GET["pic"];
  if (substr($suffix, 0, 5) == "cache"){
    $filename = $suffix;
  }
  else {
    $filename = "/tmp/MUT" . $suffix;
  }
  header('Content-Type: image/png');
  $img = imagecreatefrompng($filename);
  imagepng($img);
  imagedestroy($img);
?>
