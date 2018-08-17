while getopts ":f:g" opt; do
  case $opt in
    f)
      if ! [[ $OPTARG =~ \.fastq$ ]]
      then
        echo 'Not a fastq file'
        exit 1
      elif ! [ -f $OPTARG ]
      then
        echo 'File does not exist'
        exit 1
      fi
      ;;
    g)
      echo "-g was triggered" >&2
      ;;
    r)

      ;;
    h)
      echo "I'm helping" >&2
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done