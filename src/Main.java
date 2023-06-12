import binary.BinaryMain;
import decimal.DecimalMain;

public class Main {
    public static void main(String[] args) {

        BinaryMain binaryMain = new BinaryMain();
        DecimalMain decimalMain = new DecimalMain();
        boolean isBinary = false;
        if(isBinary){
            binaryMain.binary();
        }else{
            decimalMain.decimal();
        }
    }
}