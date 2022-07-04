function model() {
    n0 <- array();
    n1 <- array();
    n2 <- array(array(), array(), array());
    n3 <- array(array());
    n4 <- array();
    n5 <- array();
    n6 <- array(0);
    n7 <- array(0);
    n8 <- array(4, 4, 4);
    sequence_unassigned <- list(3);
    sequence_vehicle_0 <- list(3);
    timeLeavingTheWarehousevehicle_0 <- int(0, 1073741824);
    n9 <- timeLeavingTheWarehousevehicle_0 == 0;
    n10 <- array(array(0, 0, 0), array(0, 0, 0), array(0, 0, 0));
    n11 <- array(0, 0, 0);
    n12 <- array(1, 1, 1);
    n13 <- array(array(1073741824), array(1073741824), array(1073741824));
    n14 <- array(array(0), array(0), array(0));
    n15 <- array(0, 0, 0);
    n16 <- array(1, 2, 3);
    n17 <- array(0, 0, 0);
    n18 <- at(n15, at(sequence_vehicle_0, 0));
    n19 <- at(n14, at(sequence_vehicle_0, 0), 0) - at(n17, at(sequence_vehicle_0, 0));
    n20 <- int(0, 1073741824);
    n21 <- array(n20);
    n22 <- count(sequence_vehicle_0);
    n23 <- iif(n22 > 0, min(at(n21, 0), max(timeLeavingTheWarehousevehicle_0, max(0, n19 - n18))), 0);
    n24 <- array(array(1073741824), array(1073741824), array(1073741824));
    n25 <- bool();
    n26 <- bool();
    n27 <- bool();
    n28 <- array(n27, n26, n25);
    endTime_vehicle_0 <- array(0..n22 - 1, n29);
    beginTime_vehicle_0 <- array(0..n22 - 1, n30);
    arrivalTimevehicle_0 <- array(0..n22 - 1, n31);
    waitingTimevehicle_0 <- sum(0..n22 - 1, n32);
    n33 <- array(0..n22 - 1, n34);
    n35 <- sum(timeLeavingTheWarehousevehicle_0);
    n36 <- count(sequence_unassigned);
    vehicleUsed_0 <- n22 > 0;
    n37 <- vehicleUsed_0 + (n36 > 0);
    n38 <- n22 - 1;
    n39 <- array(0, 0, 0);
    routeDuration0 <- iif(n22 > 0, (at(endTime_vehicle_0, n22 - 1) + at(n39, at(sequence_vehicle_0, n38))) - n23, 0);
    totalDuration <- sum(routeDuration0);
    n40 <- n22 - 1;
    distanceToWarehouses0 <- array(3, 3, 3);
    distanceFromWarehouses0 <- array(3, 3, 3);
    routeDistance0 <- sum(1..n22 - 1, n41) + iif(n22 > 0, at(distanceFromWarehouses0, at(sequence_vehicle_0, 0)) + at(distanceToWarehouses0, at(sequence_vehicle_0, n40)), 0);
    totalDistance <- sum(routeDistance0);
    totalFixedCost <- sum(iif(vehicleUsed_0, 0, 0));
    n42 <- at(sequence_vehicle_0, n22 - 1);
    n43 <- iif(0, sum(0..0 - 1, n44), 0);
    n45 <- at(endTime_vehicle_0, n22 - 1);
    n46 <- max(iif(n22 > 0, (n45 + n43) + at(n39, n42), 0) - 2147483647, 0);
    n47 <- sum(0..n22 - 1, n34);
    n48 <- n22 > 0;
    n49 <- at(sequence_vehicle_0, n22 - 1);
    n50 <- at(distanceFromWarehouses0, at(sequence_vehicle_0, 0));
    routeDistanceCost0 <- iif(n22 > 0, ((routeDistance0 - iif(0, n50, 0)) - iif(0, at(distanceToWarehouses0, n49), 0)) * 0.0, 0);
    totalDistanceCost <- sum(routeDistanceCost0);
    n51 <- at(sequence_vehicle_0, n22 - 1);
    n52 <- at(n17, at(sequence_vehicle_0, 0));
    routeDurationCost0 <- iif(n22 > 0, ((routeDuration0 - iif(0, n52, 0)) - iif(0, at(n39, n51), 0)) * 1.0, 0);
    totalDurationCost <- sum(routeDurationCost0);
    constraint partition(sequence_vehicle_0, sequence_unassigned);
    constraint sum(sum(0..n22 - 1, n53)) == 0;
    minimize (((sum(0..n36 - 1, n54) + totalDurationCost) + totalDistanceCost) + max(0, sum(iif(n48, n47, 0) + n46 * 0.0))) + totalFixedCost;
}

function n41(arg0) {
    n55 <- arg0;
    n56 <- at(array(array(0, 3, 3), array(3, 0, 3), array(3, 3, 0)), at(sequence_vehicle_0, n55 - 1), at(sequence_vehicle_0, n55));
    return n56;
}

function n29(arg0, arg1) {
    n57 <- arg0; n58 <- arg1;
    n59 <- max(iif(n57 == 0, (n23 + at(n17, at(sequence_vehicle_0, n57))) + at(n15, at(sequence_vehicle_0, n57)), (n58 + at(n10, at(sequence_vehicle_0, n57 - 1), at(sequence_vehicle_0, n57))) + iif(at(n16, at(sequence_vehicle_0, n57)) == at(n16, at(sequence_vehicle_0, n57 - 1)), 0, at(n15, at(sequence_vehicle_0, n57)))), iif(at(n12, at(sequence_vehicle_0, n57)) == 0, 0, min(0..at(n12, at(sequence_vehicle_0, n57)) - 1, n60))) + at(n11, at(sequence_vehicle_0, n57));
    return n59;
}

function n60(arg0) {
    n61 <- arg0;
    n62 <- iif(iif(at(n28, at(sequence_vehicle_0, n57)) == 1, at(n24, at(sequence_vehicle_0, n57), n61), at(n13, at(sequence_vehicle_0, n57), n61)) >= iif(n57 == 0, (n23 + at(n17, at(sequence_vehicle_0, n57))) + at(n15, at(sequence_vehicle_0, n57)), (n58 + at(n10, at(sequence_vehicle_0, n57 - 1), at(sequence_vehicle_0, n57))) + iif(at(n16, at(sequence_vehicle_0, n57)) == at(n16, at(sequence_vehicle_0, n57 - 1)), 0, at(n15, at(sequence_vehicle_0, n57)))), at(n14, at(sequence_vehicle_0, n57), n61), at(n13, at(sequence_vehicle_0, n57), at(n12, at(sequence_vehicle_0, n57)) - 1));
    return n62;
}

function n30(arg0) {
    n63 <- arg0;
    n64 <- at(endTime_vehicle_0, n63) - at(n11, at(sequence_vehicle_0, n63));
    return n64;
}

function n31(arg0) {
    n65 <- arg0;
    n66 <- at(beginTime_vehicle_0, n65) - iif(n65 < 1, at(n15, at(sequence_vehicle_0, n65)), iif(at(n16, at(sequence_vehicle_0, n65)) == at(n16, at(sequence_vehicle_0, n65 - 1)), 0, at(n15, at(sequence_vehicle_0, n65))));
    return n66;
}

function n32(arg0) {
    n67 <- arg0;
    n68 <- iif(n67 == 0, (at(arrivalTimevehicle_0, n67) - at(n17, at(sequence_vehicle_0, n67))) - n23, max(0, (at(arrivalTimevehicle_0, n67) - at(endTime_vehicle_0, n67 - 1)) - at(n10, at(sequence_vehicle_0, n67 - 1), at(sequence_vehicle_0, n67))));
    return n68;
}

function n44(arg0) {
    n69 <- arg0;
    n70 <- iif((at(endTime_vehicle_0, n22 - 1) <= at(n21, n69)) && (at(array(2147483647), 0) >= at(n21, n69)), at(array(0), n69), 0);
    return n70;
}

function n53(arg0) {
    n71 <- arg0;
    n72 <- max(0, (at(endTime_vehicle_0, n71) - at(n11, at(sequence_vehicle_0, n71))) - at(n13, at(sequence_vehicle_0, n71), at(n12, at(sequence_vehicle_0, n71)) - 1));
    return n72;
}

function n34(arg0) {
    n73 <- arg0;
    n74 <- at(array(0.0, 0.0, 0.0), at(sequence_vehicle_0, n73)) * max(0, at(beginTime_vehicle_0, n73) - min(0..at(n12, at(sequence_vehicle_0, n73)) - 1, n75));
    return n74;
}

function n75(arg0) {
    n76 <- arg0;
    n77 <- iif(((at(n13, at(sequence_vehicle_0, n73), n76) >= at(beginTime_vehicle_0, n73)) && (at(n14, at(sequence_vehicle_0, n73), n76) <= at(beginTime_vehicle_0, n73))) && (at(n28, at(sequence_vehicle_0, n73)) == 0), at(n24, at(sequence_vehicle_0, n73), n76), at(n24, at(sequence_vehicle_0, n73), at(n12, at(sequence_vehicle_0, n73)) - 1));
    return n77;
}

function n54(arg0) {
    n78 <- arg0;
    n79 <- at(array(1000000.0, 1000000.0, 1000000.0), at(sequence_unassigned, n78));
    return n79;
}
